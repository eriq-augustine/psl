/*
 * This file is part of the PSL software.
 * Copyright 2011-2015 University of Maryland
 * Copyright 2013-2019 The Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.linqs.psl.database.rdbms;

import org.linqs.psl.config.Config;
import org.linqs.psl.database.DatabaseQuery;
import org.linqs.psl.model.atom.Atom;
import org.linqs.psl.model.formula.Conjunction;
import org.linqs.psl.model.formula.Formula;
import org.linqs.psl.model.predicate.ExternalFunctionalPredicate;
import org.linqs.psl.model.predicate.Predicate;
import org.linqs.psl.model.predicate.GroundingOnlyPredicate;
import org.linqs.psl.model.predicate.StandardPredicate;
import org.linqs.psl.model.term.Term;
import org.linqs.psl.model.term.Variable;
import org.linqs.psl.util.BitUtils;
import org.linqs.psl.util.IteratorUtils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Queue;

/**
 * Try to rewrite a grounding query to make it more efficient.
 * Of course, we would like the query to be faster (more simple) and return less results.
 * However we can't always do that, so we will tolerate more results in exchange for a faster query
 * (since trivial ground rule checking is fast).
 * Note that this class will make heavy use of referential equality.
 */
public class QueryRewriter {
    private static final Logger log = LoggerFactory.getLogger(QueryRewriter.class);

    public static final String CONFIG_PREFIX = "queryrewriter";

    /**
     * How much we allow the query cost (number of rows) to for new plans.
     */
    public static final String ALLOWED_TOTAL_INCREASE_KEY = CONFIG_PREFIX + ".allowedtotalcostincrease";
    public static final double ALLOWED_TOTAL_INCREASE_DEFAULT = 2.0;

    /**
     * How much we allow the query cost (number of rows) to increase at each step.
     */
    public static final String ALLOWED_STEP_INCREASE_KEY = CONFIG_PREFIX + ".allowedstepcostincrease";
    public static final double ALLOWED_STEP_INCREASE_DEFAULT = 1.5;

    /**
     * The different methods for estimating the join cost.
     */
    public static enum CostEstimator { SIZE, SELECTIVITY, HISTOGRAM }

    /**
     * Whether we should use histograms or column selectivity to estimate the join size.
     */
    public static final String COST_ESTIMATOR_KEY = CONFIG_PREFIX + ".costestimator";
    public static final String COST_ESTIMATOR_DEFAULT = CostEstimator.HISTOGRAM.toString();

    /**
     * The search method to use (RewriteFringe descendent).
     */
    public static final String SEARCH_TYPE_KEY = CONFIG_PREFIX + ".searchtype";
    public static final String SEARCH_TYPE_DEFAULT = BoundedRewriteFringe.class.getName();

    private RewriteFringe fringe;

    private double optimisticCostWeight;
    private double pessimisticCostWeight;
    private double optimisticRowWeight;
    private double pessimisticRowWeight;

    private double allowedTotalCostIncrease;
    private double allowedStepCostIncrease;

    private CostEstimator costEstimator;

    public QueryRewriter() {
        // TODO(eriq): Config
        optimisticCostWeight = 0.018;
        optimisticRowWeight = 0.0015;
        pessimisticCostWeight = 0.020;
        pessimisticRowWeight = 0.0020;

        allowedTotalCostIncrease = Config.getDouble(ALLOWED_TOTAL_INCREASE_KEY, ALLOWED_TOTAL_INCREASE_DEFAULT);
        allowedStepCostIncrease = Config.getDouble(ALLOWED_STEP_INCREASE_KEY, ALLOWED_STEP_INCREASE_DEFAULT);
        costEstimator = CostEstimator.valueOf(Config.getString(COST_ESTIMATOR_KEY, COST_ESTIMATOR_DEFAULT).toUpperCase());

        fringe = (RewriteFringe)Config.getNewObject(SEARCH_TYPE_KEY, SEARCH_TYPE_DEFAULT);
    }

    /**
     * Rewrite the query to minimize the execution time while trading off query size.
     */
    public Formula rewrite(Formula baseFormula, RDBMSDatabase database) {
        // Once validated, we know that the formula is a conjunction or single atom.
        DatabaseQuery.validate(baseFormula);

        // Shortcut for priors (single atoms).
        if (baseFormula instanceof Atom) {
            return baseFormula;
        }

        Set<Atom> atomBuffer = baseFormula.getAtoms(new HashSet<Atom>());
        Set<Atom> passthrough = filterBaseAtoms(atomBuffer);
        Map<Variable, Set<Atom>> variableUsageMapping = getAllUsedVariables(atomBuffer, null);

        // Get the atoms in a consistent ordering.
        List<Atom> atoms = new ArrayList<Atom>(atomBuffer);
        Collections.sort(atoms, new Comparator<Atom>() {
            public int compare(Atom a, Atom b) {
                return a.toString().compareTo(b.toString());
            }
        });

        Set<Long> seenNodes = new HashSet<Long>();
        fringe.clear();

        // A bitset representing the different atoms involved in the rewrite.
        // We will typically just pull a node's long bitset into this buffer.
        boolean[] atomBits = new boolean[atoms.size()];
        for (int i = 0; i < atomBits.length; i++) {
            atomBits[i] = true;
        }

        // Start with all atoms.
        RewriteNode baseNode = createRewriteNode(atomBits, atoms, passthrough, atomBuffer, variableUsageMapping, database);
        log.trace("Root node: " + baseNode);

        fringe.push(baseNode);
        seenNodes.add(new Long(BitUtils.toBitSet(atomBits)));

        while (fringe.size() > 0) {
            RewriteNode node = fringe.pop();

            log.trace("Expanding node: " + node);

            // Expand the node by trying to remove each atom (in-turn).
            BitUtils.toBits(node.atomsBitSet, atomBits);

            for (int i = 0; i < atoms.size(); i++) {
                // Skip if this atom was dropped somewhere else.
                if (!atomBits[i]) {
                    continue;
                }

                // Flip the atom.
                atomBits[i] = false;

                // Skip nodes we have seen before.
                Long bitId = new Long(BitUtils.toBitSet(atomBits));
                if (!seenNodes.contains(bitId)) {
                    seenNodes.add(bitId);

                    RewriteNode child = createRewriteNode(atomBits, atoms, passthrough, atomBuffer, variableUsageMapping, database);
                    log.trace("Found child: " + child);
                    fringe.push(child);
                }

                // Unflip the atom (so we can flip the next one).
                atomBits[i] = true;
            }
        }

        RewriteNode bestNode = fringe.getBestNode();
        log.debug("Computed cost-based query rewrite for [{}]({}): [{}]({}).",
                baseNode.formula, baseNode.optimisticCost, bestNode.formula, bestNode.optimisticCost);

        return bestNode.formula;
    }

    /**
     * Given the active atoms, perform the rewrite and estimate the cost.
     * The cost is non-trivial, since the query estimate is made.
     * @return The rewrite node, or null if the rewrite is invalid.
     */
    private RewriteNode createRewriteNode(boolean[] atomBits, List<Atom> atoms,
            Set<Atom> passthrough, Set<Atom> buffer, Map<Variable, Set<Atom>> variableUsageMapping,
            RDBMSDatabase database) {
        Formula formula = constructFormula(atomBits, atoms, passthrough, buffer);
        int activeAtoms = buffer.size() - passthrough.size();

        // Make sure that all variables are covered.
        for (Map.Entry<Variable, Set<Atom>> entry : variableUsageMapping.entrySet()) {
            boolean hasVariable = false;
            for (Atom atomWithVariable : entry.getValue()) {
                if (buffer.contains(atomWithVariable)) {
                    hasVariable = true;
                    break;
                }
            }

            if (!hasVariable) {
                return null;
            }
        }

        RDBMSDataStore dataStore = (RDBMSDataStore)database.getDataStore();
        String sql = Formula2SQL.getQuery(formula, database, false);
        RDBMSDataStore.ExplainResult result = dataStore.explain(sql);

        double optimisticCost = result.totalCost * optimisticCostWeight + result.rows * optimisticRowWeight * atomBits.length;
        double pessimisticCost = result.totalCost * pessimisticCostWeight + result.rows * pessimisticRowWeight * atomBits.length;

        return new RewriteNode(BitUtils.toBitSet(atomBits), formula, optimisticCost, pessimisticCost);
    }

    /**
     * Construct a formula from the active atoms.
     * |buffer| will be cleared at the beginning of this method,
     * but left filled with all the ative atoms (including passthrough) on return.
     */
    private Formula constructFormula(boolean[] atomBits, List<Atom> atoms, Set<Atom> passthrough, Set<Atom> buffer) {
        buffer.clear();
        buffer.addAll(passthrough);

        for (int i = 0; i < atomBits.length; i++) {
            if (atomBits[i]) {
                buffer.add(atoms.get(i));
            }
        }

        Formula formula = null;
        if (buffer.size() == 1) {
            formula = buffer.iterator().next();
        } else {
            formula = new Conjunction(buffer.toArray(new Formula[0]));
        }

        return formula;
    }

    /**
     * Rewrite the query to minimize the execution time acording to join estimators while trading off query size.
     */
    public Formula rewriteWithJoinEstimation(Formula baseFormula, RDBMSDataStore dataStore) {
        // Once validated, we know that the formula is a conjunction or single atom.
        DatabaseQuery.validate(baseFormula);

        // Shortcut for priors (single atoms).
        if (baseFormula instanceof Atom) {
            return baseFormula;
        }

        Set<Atom> usedAtoms = baseFormula.getAtoms(new HashSet<Atom>());
        Set<Atom> passthrough = filterBaseAtoms(usedAtoms);

        Map<Predicate, TableStats> tableStats = fetchTableStats(usedAtoms, dataStore);

        double baseCost = estimateQuerySize(costEstimator, usedAtoms, null, tableStats, dataStore);
        double currentCost = baseCost;

        log.trace("Starting cost: " + baseCost);

        // TODO(eriq): A DP approach should be better than a greedy one.

        while (true) {
            // The estimated query cost after removing the target atom.
            double bestCost = -1;
            Atom bestAtom = null;

            for (Atom atom : usedAtoms) {
                if (canRemove(atom, usedAtoms)) {
                    double cost = estimateQuerySize(costEstimator, usedAtoms, atom, tableStats, dataStore);
                    if (cost < 0.0) {
                        log.trace("Planned Cost for (" + usedAtoms + " - " + atom + "): MAX");
                        continue;
                    }

                    log.trace("Planned Cost for (" + usedAtoms + " - " + atom + "): " + cost);

                    if (bestAtom == null || cost < bestCost) {
                        bestAtom = atom;
                        bestCost = cost;
                    }
                }
            }

            if (bestAtom == null) {
                break;
            }

            // We couldn't find any viable plans.
            if (bestAtom == null) {
                break;
            }

            // We expect the cost to go up, but will cut it off at some point.
            if (bestCost > (baseCost * allowedTotalCostIncrease) || bestCost > (currentCost * allowedStepCostIncrease)) {
                break;
            }

            usedAtoms.remove(bestAtom);
            currentCost = bestCost;

            log.trace("Choose plan for iteration: " + usedAtoms + ": " + bestCost);
        }

        usedAtoms.addAll(passthrough);

        Formula query = null;
        if (usedAtoms.size() == 1) {
            query = usedAtoms.iterator().next();
        } else {
            query = new Conjunction(usedAtoms.toArray(new Formula[0]));
        }

        log.debug("Computed cost-based query rewrite for [{}]({}): [{}]({}).", baseFormula, baseCost, query, currentCost);
        return query;
    }

    /**
     * Get all the valid candidate queries for a base query.
     * Rules with more than 64 standard atoms will throw.
     */
    public List<Formula> allCandidates(Formula baseFormula) {
        // Once validated, we know that the formula is a conjunction or single atom.
        DatabaseQuery.validate(baseFormula);

        Set<Atom> baseAtoms = baseFormula.getAtoms(new HashSet<Atom>());
        Set<Atom> passthroughAtoms = filterBaseAtoms(baseAtoms);

        List<Atom> atoms = new ArrayList<Atom>(baseAtoms);
        List<Formula> candidates = new ArrayList<Formula>((int)Math.pow(2, atoms.size()));

        Set<Variable> allVariables = new HashSet<Variable>();
        Map<Atom, Set<Variable>> atomVariables = new HashMap<Atom, Set<Variable>>(atoms.size());
        for (Atom atom : atoms) {
            Set<Variable> variables = atom.getVariables();

            allVariables.addAll(variables);
            atomVariables.put(atom, variables);
        }

        Set<Variable> usedVariables = new HashSet<Variable>(allVariables.size());
        List<Atom> candidate = new ArrayList<Atom>(atoms.size());

        for (boolean[] activeAtoms : IteratorUtils.powerset(atoms.size())) {
            usedVariables.clear();
            candidate.clear();

            for (int i = 0; i < atoms.size(); i++) {
                if (activeAtoms[i]) {
                    usedVariables.addAll(atomVariables.get(atoms.get(i)));
                    candidate.add(atoms.get(i));
                }
            }

            // If all variables were not covered.
            if (usedVariables.size() != allVariables.size()) {
                continue;
            }

            candidate.addAll(passthroughAtoms);

            Formula query = null;
            if (candidate.size() == 1) {
                query = candidate.get(0);
            } else {
                query = new Conjunction(candidate.toArray(new Formula[0]));
            }

            candidates.add(query);
        }

        return candidates;
    }

    private double estimateQuerySize(CostEstimator costEstimator, Set<Atom> atoms, Atom ignore, Map<Predicate, TableStats> tableStats, RDBMSDataStore dataStore) {
        if (costEstimator == CostEstimator.HISTOGRAM) {
            return estimateQuerySizeWithHistorgram(atoms, ignore, tableStats, dataStore);
        } else if (costEstimator == CostEstimator.SELECTIVITY) {
            return estimateQuerySizeWithSelectivity(atoms, ignore, tableStats, dataStore);
        } else if (costEstimator == CostEstimator.SIZE) {
            return estimateQuerySizeWithSize(atoms, ignore, tableStats, dataStore);
        }

        throw new IllegalStateException("Unknown CostEstimator value: " + costEstimator);
    }

    /**
     * Estimate the cost of the query (conjunctive query over the given atoms) using histogram stats.
     * Based off of: http://consystlab.unl.edu/Documents/StudentReports/Working-Note2-2008.pdf
     * @param ignore if not null, then do not include it in the cost computation.
     * @return negative if we exceed max size and should just discard this plan.
     */
    private double estimateQuerySizeWithHistorgram(Set<Atom> atoms, Atom ignore, Map<Predicate, TableStats> tableStats, RDBMSDataStore dataStore) {
        double cost = 1.0;

        // Start with the product of the joins.
        for (Atom atom : atoms) {
            if (atom == ignore) {
                continue;
            }

            cost *= tableStats.get(atom.getPredicate()).getCount();
        }

        // Now compute for each variable.
        for (Map.Entry<Variable, Set<Atom>> entry : getAllUsedVariables(atoms, null).entrySet()) {
            Variable variable = entry.getKey();
            Set<Atom> involvedAtoms = entry.getValue();

            // Only worry about join columns.
            if (involvedAtoms.size() <= 1) {
                continue;
            }

            // The histogram that we will chain together using all the involved atoms.
            SelectivityHistogram joinHistogram = null;
            double crossProductSize = 1.0;

            for (Atom atom : involvedAtoms) {
                if (atom == ignore) {
                    continue;
                }

                crossProductSize *= tableStats.get(atom.getPredicate()).getCount();

                String columnName = getColumnName(dataStore, atom, variable);
                SelectivityHistogram histogram = tableStats.get(atom.getPredicate()).getHistogram(columnName);

                if (joinHistogram == null) {
                    joinHistogram = histogram;
                } else {
                    @SuppressWarnings("unchecked")
                    SelectivityHistogram ignoreUnchecked = joinHistogram.join(histogram);
                    joinHistogram = ignoreUnchecked;
                }
            }

            long joinSize = joinHistogram.size();
            if (joinSize == Long.MAX_VALUE) {
                return -1.0;
            }

            cost *= (joinHistogram.size() / crossProductSize);
        }

        return cost;
    }

    /**
     * Estimate the cost of the query (conjunctive query over the given atoms).
     * Based off of: http://users.csc.calpoly.edu/~dekhtyar/468-Spring2016/lectures/lec17.468.pdf
     * @param ignore if not null, then do not include it in the cost computation.
     */
    private double estimateQuerySizeWithSelectivity(Set<Atom> atoms, Atom ignore, Map<Predicate, TableStats> tableStats, RDBMSDataStore dataStore) {
        double cost = 1.0;

        // Start with the product of the joins.
        for (Atom atom : atoms) {
            if (atom == ignore) {
                continue;
            }

            cost *= tableStats.get(atom.getPredicate()).getCount();
        }

        // Now take out the join factor for each join attribute.
        for (Variable variable : getAllUsedVariables(atoms, null).keySet()) {
            int cardinalitiesCount = 0;
            int minCardinality = 0;

            for (Atom atom : atoms) {
                if (atom == ignore) {
                    continue;
                }

                if (!atom.getVariables().contains(variable)) {
                    continue;
                }

                String columnName = getColumnName(dataStore, atom, variable);
                int cardinality = tableStats.get(atom.getPredicate()).getCardinality(columnName);

                if (cardinalitiesCount == 0 || cardinality < minCardinality) {
                    minCardinality = cardinality;
                }
                cardinalitiesCount += 1;
            }

            // Only worry about join columns.
            if (cardinalitiesCount <= 1) {
                continue;
            }

            cost /= (double)minCardinality;
        }

        return cost;
    }

    /**
     * Estimate the cost of the query (conjunctive query over the given atoms) using just the size of the involved tables.
     * @param ignore if not null, then do not include it in the cost computation.
     */
    private double estimateQuerySizeWithSize(Set<Atom> atoms, Atom ignore, Map<Predicate, TableStats> tableStats, RDBMSDataStore dataStore) {
        double cost = 1.0;

        // Just use the product of the joins.
        for (Atom atom : atoms) {
            if (atom == ignore) {
                continue;
            }

            cost *= tableStats.get(atom.getPredicate()).getCount();
        }

        return cost;
    }

    private String getColumnName(RDBMSDataStore dataStore, Atom atom, Variable variable) {
        int columnIndex = -1;
        Term[] args = atom.getArguments();
        for (int i = 0; i < args.length; i++) {
            if (variable.equals(args[i])) {
                columnIndex = i;
                break;
            }
        }

        if (columnIndex == -1) {
            throw new java.util.NoSuchElementException(String.format("Could not find column name for variable %s in atom %s.", variable, atom));
        }

        return dataStore.getPredicateInfo(atom.getPredicate()).argumentColumns().get(columnIndex);
    }

    /**
     * Is it safe to remove the given atom from the given query?
     */
    private boolean canRemove(Atom atom, Set<Atom> usedAtoms) {
        // Make sure that we do not remove any variables (ie there is at least one other atom that uses the variable)..
        Set<Variable> remainingVariables = atom.getVariables();
        remainingVariables.removeAll(getAllUsedVariables(usedAtoms, atom).keySet());

        return remainingVariables.size() == 0;
    }

    private Map<Variable, Set<Atom>> getAllUsedVariables(Set<Atom> atoms, Atom ignore) {
        Map<Variable, Set<Atom>> variables = new HashMap<Variable, Set<Atom>>();

        for (Atom atom : atoms) {
            if (atom == ignore) {
                continue;
            }

            for (Variable variable : atom.getVariables()) {
                if (!variables.containsKey(variable)) {
                    variables.put(variable, new HashSet<Atom>());
                }

                variables.get(variable).add(atom);
            }
        }

        return variables;
    }

    private Map<Predicate, TableStats> fetchTableStats(Set<Atom> usedAtoms, RDBMSDataStore dataStore) {
        Set<Predicate> predicates = new HashSet<Predicate>();
        for (Atom atom : usedAtoms) {
            predicates.add(atom.getPredicate());
        }

        Map<Predicate, TableStats> stats = new HashMap<Predicate, TableStats>();
        for (Predicate predicate : predicates) {
            stats.put(predicate, dataStore.getPredicateInfo(predicate).getTableStats(dataStore.getDriver()));
        }

        return stats;
    }

    /**
     * Filter the initial set of atoms.
     * Remove external functional prediates and pass through grounding only predicates.
    */
    private Set<Atom> filterBaseAtoms(Set<Atom> atoms) {
        Set<Atom> passthrough = new HashSet<Atom>();

        Set<Atom> removeAtoms = new HashSet<Atom>();
        for (Atom atom : atoms) {
            if (atom.getPredicate() instanceof ExternalFunctionalPredicate) {
                // Skip. These are handled at instantiation time.
                removeAtoms.add(atom);
            } else if (atom.getPredicate() instanceof GroundingOnlyPredicate) {
                // Passthrough.
                removeAtoms.add(atom);
                passthrough.add(atom);
            } else if (!(atom.getPredicate() instanceof StandardPredicate)) {
                throw new IllegalStateException("Unknown predicate type: " + atom.getPredicate().getClass().getName());
            }
        }

        atoms.removeAll(removeAtoms);
        return passthrough;
    }

    private static class RewriteNode {
        public long atomsBitSet;
        public Formula formula;
        public double optimisticCost;
        public double pessimisticCost;

        public RewriteNode(long atomsBitSet, Formula formula, double optimisticCost, double pessimisticCost) {
            this.atomsBitSet = atomsBitSet;
            this.formula = formula;
            this.optimisticCost = optimisticCost;
            this.pessimisticCost = pessimisticCost;
        }

        public String toString() {
            return String.format("{Atom Bits: %d, Optimistic: %f, Pessimistic: %f, Formula: %s}",
                    atomsBitSet, optimisticCost, pessimisticCost, formula);
        }
    }

    /**
     * Defines the strategy for exploring the rewrite space.
     */
    public abstract static class RewriteFringe {
        protected RewriteNode bestNode;
        protected Collection<RewriteNode> fringe;

        protected RewriteFringe(Collection<RewriteNode> fringe) {
            bestNode = null;
            this.fringe = fringe;
        }

        public int size() {
            return fringe.size();
        }

        public void clear() {
            bestNode = null;
            fringe.clear();
        }

        /**
         * Add a node to the fringe.
         * The node will be check if it is the best on add because the search through
         * the rewrite space may be time/cost bounded,
         *  and we want to check for the best node as soon as we pay the estimation cost (SQL EXPLAIN).
         */
        public void push(RewriteNode node) {
            // Skip invalid rewrites.
            if (node == null) {
                return;
            }

            if (bestNode == null || node.optimisticCost < bestNode.optimisticCost) {
                bestNode = node;
            }

            if (pushInternal(node)) {
                log.trace("   Accepted node.");
            } else {
                log.trace("   Rejected node.");
            }
        }

        public RewriteNode getBestNode() {
            return bestNode;
        }

        /**
         * Actually add the node to the fringe.
         * @return true is the node was accepted, false otherwise.
         */
        protected abstract boolean pushInternal(RewriteNode node);
        public abstract RewriteNode pop();
    }

    public static class DFSRewriteFringe extends RewriteFringe {
        public DFSRewriteFringe() {
            super(new LinkedList<RewriteNode>());
        }

        protected boolean pushInternal(RewriteNode node) {
            ((LinkedList<RewriteNode>)fringe).addFirst(node);
            return true;
        }

        public RewriteNode pop() {
            return ((LinkedList<RewriteNode>)fringe).removeFirst();
        }
    }

    public static class BFSRewriteFringe extends RewriteFringe {
        public BFSRewriteFringe() {
            super(new LinkedList<RewriteNode>());
        }

        protected boolean pushInternal(RewriteNode node) {
            ((LinkedList<RewriteNode>)fringe).addLast(node);
            return true;
        }

        public RewriteNode pop() {
            return ((LinkedList<RewriteNode>)fringe).removeFirst();
        }
    }

    public static class UCSRewriteFringe extends RewriteFringe {
        public UCSRewriteFringe() {
            // Use a priority queue (with PriorityQueue's default size of 11.
            super(new PriorityQueue<RewriteNode>(11, new Comparator<RewriteNode>() {
                public int compare(RewriteNode a, RewriteNode b) {
                    if (a.optimisticCost < b.optimisticCost) {
                        return -1;
                    }

                    if (a.optimisticCost > b.optimisticCost) {
                        return 1;
                    }

                    return 0;
                }
            }));
        }

        protected boolean pushInternal(RewriteNode node) {
            ((PriorityQueue<RewriteNode>)fringe).add(node);
            return true;
        }

        public RewriteNode pop() {
            return ((PriorityQueue<RewriteNode>)fringe).poll();
        }
    }

    public static class BoundedRewriteFringe extends UCSRewriteFringe {
        @Override
        protected boolean pushInternal(RewriteNode node) {
            if (node.optimisticCost > bestNode.pessimisticCost) {
                return false;
            }

            return super.pushInternal(node);
        }
    }
}
