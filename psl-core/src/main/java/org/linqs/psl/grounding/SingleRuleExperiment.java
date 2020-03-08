/*
 * This file is part of the PSL software.
 * Copyright 2011-2015 University of Maryland
 * Copyright 2013-2020 The Regents of the University of California
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
package org.linqs.psl.grounding;

import org.linqs.psl.config.Config;
import org.linqs.psl.database.DataStore;
import org.linqs.psl.database.atom.AtomManager;
import org.linqs.psl.database.rdbms.Formula2SQL;
import org.linqs.psl.database.rdbms.RDBMSDataStore;
import org.linqs.psl.database.rdbms.RDBMSDatabase;
import org.linqs.psl.grounding.rewrite.QueryRewriter;
import org.linqs.psl.model.atom.Atom;
import org.linqs.psl.model.formula.Formula;
import org.linqs.psl.model.predicate.StandardPredicate;
import org.linqs.psl.model.rule.Rule;
import org.linqs.psl.model.rule.logical.AbstractLogicalRule;
import org.linqs.psl.model.term.ConstantType;
import org.linqs.psl.model.term.Term;
import org.linqs.psl.model.term.Variable;
import org.linqs.psl.util.ListUtils;
import org.linqs.psl.util.Parallel;
import org.linqs.psl.util.StringUtils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Run a grounding experiment where we run a single rule at a time and collect performance stats.
 */
public class SingleRuleExperiment {
    private static final Logger log = LoggerFactory.getLogger(SingleRuleExperiment.class);

    public static final String CONFIG_PREFIX = "grounding.experiment";

    public static final String EXPERIMENT_KEY = CONFIG_PREFIX;
    public static final boolean EXPERIMENT_DEFAULT = false;

    /**
     * The specific rules and rewrites to run.
     * "rule:rewrite;...".
     * QUERY_VIEW_ALL_VALUE if you want to get all the rewrites printed.
     * RULE_QUERIES_DEFAULT if you want to get the number of rules.
     */
    public static final String RULE_QUERIES_KEY = CONFIG_PREFIX + ".rulequeries";
    public static final String RULE_QUERIES_DEFAULT = "";

    /**
     * For each rule, run a full estimation (descend the entire search tree).
     */
    public static final String FULL_ESTIMATE_KEY = CONFIG_PREFIX + ".fullestimate";
    public static final boolean FULL_ESTIMATE_DEFAULT = false;

    public static final String RULE_DELIM = "_";
    public static final String QUERY_DELIM = ":";

    public static final int QUERY_VIEW_ALL_VALUE = -1;

    // Static only.
    private SingleRuleExperiment() {}

    /**
     * Run rewites of specific rules (chosen by RULE_QUERIES_KEY).
     */
    public static void run(List<Rule> rules, AtomManager atomManager, GroundRuleStore groundRuleStore) {
        String ruleQueriesString = Config.getString(RULE_QUERIES_KEY, RULE_QUERIES_DEFAULT);

        if (ruleQueriesString.equals(RULE_QUERIES_DEFAULT)) {
            log.info("Total available rules: {}", rules.size());
            return;
        }

        boolean fullEstimate = Config.getBoolean(FULL_ESTIMATE_KEY, FULL_ESTIMATE_DEFAULT);
        String[] ruleQueries = ruleQueriesString.split(RULE_DELIM);

        Parallel.initPool();

        log.info("Single rule experiment total available rules: {}", rules.size());

        if (ruleQueries.length == 0) {
            return;
        }

        DataStore dataStore = atomManager.getDatabase().getDataStore();
        QueryRewriter rewriter = new QueryRewriter();

        for (String ruleQuery : ruleQueries) {
            int[] values = StringUtils.splitInt(ruleQuery, QUERY_DELIM);
            int ruleIndex = values[0];
            int queryIndex = values[1];

            Rule rule = rules.get(ruleIndex);

            runSingleRule(rule, ruleIndex, queryIndex, atomManager, groundRuleStore, dataStore, rewriter, fullEstimate);
        }
    }

    private static void runSingleRule(Rule rule, int ruleIndex, int queryToGround,
            AtomManager atomManager, GroundRuleStore groundRuleStore,
            DataStore dataStore, QueryRewriter rewriter,
            boolean fullEstimate) {
        log.info("Single rule experiment on rule {} -- {}", ruleIndex, rule);

        List<Formula> queries = rewriter.allCandidates(rule.getRewritableGroundingFormula(atomManager));
        log.info("Found {} candidate queries.", queries.size());

        if (queryToGround != QUERY_VIEW_ALL_VALUE && (queryToGround < 0 || queryToGround > queries.size())) {
            throw new RuntimeException("Bad value for query (rewrite), got: " + queryToGround);
        }

        // Get a consistent ordering of all the queries.
        List<String> queryKeys = new ArrayList<String>();
        Map<String, Formula> queryMapping = new HashMap<String, Formula>();

        buildQueryKeys(queries, queryKeys, queryMapping);

        if (fullEstimate) {
            log.info("Estimating full search space.");
            List<String> row = new ArrayList<String>();

            row.add("Index");
            row.add("Query");
            row.add("Count");
            row.add("Cost");
            row.add("Rows");
            log.info("FullEstimate -- " + ListUtils.join("\t", row));

            for (int i = 0; i < queryKeys.size(); i++) {
                row.clear();

                String sql = Formula2SQL.getQuery(queryMapping.get(queryKeys.get(i)), (RDBMSDatabase)atomManager.getDatabase(), false);
                RDBMSDataStore.ExplainResult result = ((RDBMSDataStore)dataStore).explain(sql);

                row.add("" + i);
                row.add(queryMapping.get(queryKeys.get(i)).toString());
                row.add("" + result.totalCost);
                row.add("" + result.rows);
                log.info("FullEstimate -- " + ListUtils.join("\t", row));
            }
        }

        if (queryToGround == QUERY_VIEW_ALL_VALUE) {
            for (int i = 0; i < queryKeys.size(); i++) {
                log.info("   {} -- {}", i, queryMapping.get(queryKeys.get(i)));
            }

            return;
        }

        List<String> header = new ArrayList<String>();
        List<String> stats = new ArrayList<String>();

        Formula query = queryMapping.get(queryKeys.get(queryToGround));

        addStats(header, stats, "Rule Index", ruleIndex);
        addStats(header, stats, "Rule String", rule);
        addStats(header, stats, "Query Index", queryToGround);
        addStats(header, stats, "Query String", query);
        addStats(header, stats, "Logical?", "" + (rule instanceof AbstractLogicalRule));

        computeStaticStats(header, stats, query);

        String sql = Formula2SQL.getQuery(query, (RDBMSDatabase)atomManager.getDatabase(), false);
        RDBMSDataStore.ExplainResult explainResult = ((RDBMSDataStore)dataStore).explain(sql);

        addStats(header, stats, "Estimated Startup Cost", explainResult.startupCost);
        addStats(header, stats, "Estimated Total Cost", explainResult.totalCost);
        addStats(header, stats, "Estimated Rows", explainResult.rows);

        List<Rule> tempRules = new ArrayList<Rule>();
        tempRules.add(rule);
        Grounding.groundParallel(query, tempRules, atomManager, groundRuleStore, false);

        log.info("Query Stats Header: " + ListUtils.join("\t", header));
        log.info("Query Stats Data: " + ListUtils.join("\t", stats));
    }

    private static void buildQueryKeys(List<Formula> queries, List<String> queryKeys, Map<String, Formula> queryMapping) {
        for (Formula query : queries) {
            String key = StringUtils.sort(query.toString());

            int atomCount = 0;
            for (Atom atom : query.getAtoms(new HashSet<Atom>())) {
                if (atom.getPredicate() instanceof StandardPredicate) {
                    atomCount++;
                }
            }

            key = String.format("%03d -- %s", atomCount, key);

            queryKeys.add(key);
            queryMapping.put(key, query);
        }

        Collections.sort(queryKeys);
        Collections.reverse(queryKeys);
    }

    /**
     * Compute stats about the static attributes of a formula.
     */
    private static void computeStaticStats(List<String> header, List<String> stats, Formula query) {
        Map<Variable, Integer> atomsPerVariable = new HashMap<Variable, Integer>();
        int numAtoms = 0;
        int numVariables = 0;
        int numUniqueIntVariables = 0;

        Set<Atom> atoms = new HashSet<Atom>();
        query.getAtoms(atoms);

        // Go through every variable.
        for (Atom atom : atoms) {
            if (!(atom.getPredicate() instanceof StandardPredicate)) {
                continue;
            }

            numAtoms++;

            Term[] args = atom.getArguments();
            for (int i = 0; i < args.length; i++) {
                if (!(args[i] instanceof Variable)) {
                    continue;
                }

                Variable variable = (Variable)args[i];
                numVariables++;

                if (atom.getPredicate().getArgumentType(i) == ConstantType.UniqueIntID) {
                    numUniqueIntVariables++;
                }

                if (!atomsPerVariable.containsKey(variable)) {
                    atomsPerVariable.put(variable, new Integer(0));
                }

                atomsPerVariable.put(variable, new Integer(atomsPerVariable.get(variable).intValue() + 1));
            }
        }

        double meanAtomsPerVariable = 0.0;
        for (Integer count : atomsPerVariable.values()) {
            meanAtomsPerVariable += count;
        }
        meanAtomsPerVariable /= atomsPerVariable.size();

        addStats(header, stats, "Atom Count", numAtoms);
        addStats(header, stats, "Variable Count", numVariables);
        addStats(header, stats, "Unique Variable Count", atomsPerVariable.size());
        addStats(header, stats, "Percent Unique Int Variables", numUniqueIntVariables / (double)numVariables);
        addStats(header, stats, "Percent Other Variables", (numVariables - numUniqueIntVariables) / (double)numVariables);
        addStats(header, stats, "Mean Atoms Per Variable", meanAtomsPerVariable);
    }

    private static void addStats(List<String> headers, List<String> stats, String header, int stat) {
        addStats(headers, stats, header, "" + stat);
    }

    private static void addStats(List<String> headers, List<String> stats, String header, double stat) {
        addStats(headers, stats, header, "" + stat);
    }

    private static void addStats(List<String> headers, List<String> stats, String header, Object stat) {
        headers.add(header);
        stats.add(stat.toString());
    }
}
