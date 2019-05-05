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
package org.linqs.psl.application.util;

import org.linqs.psl.application.groundrulestore.GroundRuleStore;
import org.linqs.psl.config.Config;
import org.linqs.psl.database.DataStore;
import org.linqs.psl.database.QueryResultIterable;
import org.linqs.psl.database.atom.AtomManager;
import org.linqs.psl.database.rdbms.Formula2SQL;
import org.linqs.psl.database.rdbms.QueryRewriter;
import org.linqs.psl.database.rdbms.RDBMSDataStore;
import org.linqs.psl.database.rdbms.RDBMSDatabase;
import org.linqs.psl.model.Model;
import org.linqs.psl.model.atom.Atom;
import org.linqs.psl.model.formula.Formula;
import org.linqs.psl.model.predicate.StandardPredicate;
import org.linqs.psl.model.rule.GroundRule;
import org.linqs.psl.model.rule.Rule;
import org.linqs.psl.model.term.Constant;
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
 * Static utilities for common {@link Model}-grounding tasks.
 */
public class Grounding {
    private static final Logger log = LoggerFactory.getLogger(Grounding.class);

    public static final String CONFIG_PREFIX = "grounding";

    /**
     * Potentially rewrite the grounding queries.
     */
    public static final String REWRITE_QUERY_KEY = CONFIG_PREFIX + ".rewritequeries";
    public static final boolean REWRITE_QUERY_DEFAULT = false;

    /**
     * Whether or not queries are being rewritten, perform the grounding queries one at a time.
     */
    public static final String SERIAL_KEY = CONFIG_PREFIX + ".serial";
    public static final boolean SERIAL_DEFAULT = false;

    /**
     * Whether or not to start instantiating ground rules as query results stream in,
     * instead of waiting for them all the show up.
     * TODO(eriq): Will be true in real code.
     */
    public static final String EAGER_INSTANTIATION_KEY = CONFIG_PREFIX + ".eagerinstantiation";
    public static final boolean EAGER_INSTANTIATION_DEFAULT = false;

    public static final String EXPERIMENT_KEY = CONFIG_PREFIX + ".experiment";
    public static final boolean EXPERIMENT_DEFAULT = false;

    /**
     * The specific rules and rewrites to run.
     * "rule:rewrite;...".
     * QUERY_VIEW_ALL_VALUE value if you want to get all the rewrites printed.
     */
    public static final String EXPERIMENT_RULE_QUERIES_KEY = CONFIG_PREFIX + ".experiment.rulequeries";
    public static final String EXPERIMENT_RULE_QUERIES_DEFAULT = "";

    /**
     * For each rule, run a full estimation (descend the entire search tree).
     */
    public static final String EXPERIMENT_FULL_ESTIMATE_KEY = CONFIG_PREFIX + ".experiment.fullestimate";
    public static final boolean EXPERIMENT_FULL_ESTIMATE_DEFAULT = false;

    public static final String RULE_DELIM = "_";
    public static final String QUERY_DELIM = ":";
    public static final int QUERY_VIEW_ALL_VALUE = -1;

    public static final int PARALLEL_WORKER_QUEUE_SIZE = 1000;

    // Static only.
    private Grounding() {}

    /**
     * Ground all the given rules.
     * @return the number of ground rules generated.
     */
    public static int groundAll(Model model, AtomManager atomManager, GroundRuleStore groundRuleStore) {
        return groundAll(model.getRules(), atomManager, groundRuleStore);
    }

    /**
     * Ground all the given rules one at a time.
     * Callers should prefer groundAll() to this since it will perform a more efficient grounding.
     * @return the number of ground rules generated.
     */
    public static int groundAllSerial(List<Rule> rules, AtomManager atomManager, GroundRuleStore groundRuleStore) {
        int groundCount = 0;
        for (Rule rule : rules) {
            groundCount += rule.groundAll(atomManager, groundRuleStore);
        }

        return groundCount;
    }

    /**
     * Ground all the given rules.
     * @return the number of ground rules generated.
     */
    public static int groundAll(List<Rule> rules, AtomManager atomManager, GroundRuleStore groundRuleStore) {
        boolean rewrite = Config.getBoolean(REWRITE_QUERY_KEY, REWRITE_QUERY_DEFAULT);
        boolean serial = Config.getBoolean(SERIAL_KEY, SERIAL_DEFAULT);
        boolean eagerInstantiation = Config.getBoolean(EAGER_INSTANTIATION_KEY, EAGER_INSTANTIATION_DEFAULT);

        Map<Formula, List<Rule>> queries = new HashMap<Formula, List<Rule>>();
        List<Rule> bypassRules = new ArrayList<Rule>();

        DataStore dataStore = atomManager.getDatabase().getDataStore();
        if (rewrite && !(dataStore instanceof RDBMSDataStore)) {
            log.warn("Cannot rewrite queries with a non-RDBMS DataStore. Queries will not be rewritten.");
            rewrite = false;
        }

        QueryRewriter rewriter = null;
        if (rewrite) {
            rewriter = new QueryRewriter();
        }

        for (Rule rule : rules) {
            if (!rule.supportsIndividualGrounding()) {
                bypassRules.add(rule);
                continue;
            }

            Formula query = rule.getGroundingFormula();
            if (rewrite) {
                query = rewriter.rewrite(query, (RDBMSDatabase)atomManager.getDatabase());
            }

            if (!queries.containsKey(query)) {
                queries.put(query, new ArrayList<Rule>());
            }

            queries.get(query).add(rule);
        }

        int initialSize = groundRuleStore.size();

        // First perform all the rewritten querties.
        for (Map.Entry<Formula, List<Rule>> entry : queries.entrySet()) {
            if (!serial) {
                // If parallel, ground all the rules that match this formula at once.
                groundParallel(entry.getKey(), entry.getValue(), atomManager, groundRuleStore, eagerInstantiation);
            } else {
                // If serial, ground the rules with this formula one at a time.
                for (Rule rule : entry.getValue()) {
                    List<Rule> tempRules = new ArrayList<Rule>();
                    tempRules.add(rule);

                    groundParallel(entry.getKey(), tempRules, atomManager, groundRuleStore, eagerInstantiation);
                }
            }
        }

        // Now ground the bypassed rules.
        groundAllSerial(bypassRules, atomManager, groundRuleStore);

        return groundRuleStore.size() - initialSize;
    }

    /**
     * Experiment only.
     * Run rewites of specific rules (chosen by EXPERIMENT_RULE_QUERIES_KEY).
     */
    public static void groundingExperiment(List<Rule> rules, AtomManager atomManager, GroundRuleStore groundRuleStore) {
        String ruleQueriesString = Config.getString(EXPERIMENT_RULE_QUERIES_KEY, EXPERIMENT_RULE_QUERIES_DEFAULT);
        boolean fullEstimate = Config.getBoolean(EXPERIMENT_FULL_ESTIMATE_KEY, EXPERIMENT_FULL_ESTIMATE_DEFAULT);
        String[] ruleQueries = ruleQueriesString.split(RULE_DELIM);

        Parallel.initPool();

        log.info("Grounding experiment total available rules: {}", rules.size());

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

            singleGroundingExperiment(rule, ruleIndex, queryIndex, atomManager, groundRuleStore, dataStore, rewriter, fullEstimate);
        }
    }

    private static void singleGroundingExperiment(Rule rule, int ruleIndex, int queryToGround,
            AtomManager atomManager, GroundRuleStore groundRuleStore,
            DataStore dataStore, QueryRewriter rewriter,
            boolean fullEstimate) {
        log.info("Grounding experiment on rule {} -- {}", ruleIndex, rule);

        List<Formula> queries = rewriter.allCandidates(rule.getGroundingFormula());

        // Get a consistent ordering of all the queries.
        List<String> queryKeys = new ArrayList<String>();
        Map<String, Formula> queryMapping = new HashMap<String, Formula>();
        Map<String, Integer> atomCounts = new HashMap<String, Integer>();

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
            atomCounts.put(key, atomCount);
        }

        Collections.sort(queryKeys);
        Collections.reverse(queryKeys);

        log.info("Found {} candidate queries.", queries.size());

        if (queryToGround != QUERY_VIEW_ALL_VALUE && (queryToGround < 0 || queryToGround > queries.size())) {
            throw new RuntimeException("Bad value for query (rewrite), got: " + queryToGround);
        }

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
                row.add(atomCounts.get(queryKeys.get(i)).toString());
                row.add("" + result.totalCost);
                row.add("" + result.rows);
                log.info("FullEstimate -- " + ListUtils.join("\t", row));
            }
        }

        if (queryToGround == QUERY_VIEW_ALL_VALUE) {
            for (int i = 0; i < queryKeys.size(); i++) {
                log.info("   {} -- {}", i, queryMapping.get(queryKeys.get(i)));
            }
        } else {
            Formula query = queryMapping.get(queryKeys.get(queryToGround));
            Integer atomCount = atomCounts.get(queryKeys.get(queryToGround));

            log.info("Query {} -- Formula: {}", queryToGround, query);
            log.info("Query {} -- Atom Count: {}", queryToGround, atomCount);

            List<Rule> tempRules = new ArrayList<Rule>();
            tempRules.add(rule);

            groundParallel(query, tempRules, atomManager, groundRuleStore, false);
        }
    }

    private static int groundParallel(Formula query, List<Rule> rules, AtomManager atomManager, GroundRuleStore groundRuleStore, boolean eagerInstantiation) {
        log.debug("Grounding {} rule(s) with query: [{}].", rules.size(), query);
        for (Rule rule : rules) {
            log.trace("    " + rule);
        }

        // We will manually handle these in the grounding process.
        // We do not want to throw too early because the ground rule may turn out to be trivial in the end.
        boolean oldAccessExceptionState = atomManager.enableAccessExceptions(false);

        int initialCount = groundRuleStore.size();

        // NOTE(eriq): For experiments, don't instantiate as we receive results in.
        //  Instead, read them all into memory and then call after the query complete.
        QueryResultIterable queryResults = atomManager.executeGroundingQuery(query);
        if (eagerInstantiation) {
            Parallel.RunTimings timings = Parallel.foreach(queryResults, new GroundWorker(atomManager, groundRuleStore, queryResults.getVariableMap(), rules));
            log.trace("Got {} results from query [{}].", timings.iterations, query);
        } else {
            List<Constant[]> results = new ArrayList<Constant[]>();
            boolean first = true;

            for (Constant[] tuple : queryResults) {
                if (first) {
                    first = false;
                    log.info("First Query Response");
                }

                results.add(tuple);
            }

            log.info("Query Complete");
            log.trace("Got {} results from query [{}].", results.size(), query);

            Parallel.foreach(results, PARALLEL_WORKER_QUEUE_SIZE, new GroundWorker(atomManager, groundRuleStore, queryResults.getVariableMap(), rules));
        }

        int groundCount = groundRuleStore.size() - initialCount;

        atomManager.enableAccessExceptions(oldAccessExceptionState);

        log.debug("Generated {} ground rules with query: [{}].", groundCount, query);
        return groundCount;
    }

    private static class GroundWorker extends Parallel.Worker<Constant[]> {
        private AtomManager atomManager;
        private GroundRuleStore groundRuleStore;
        private Map<Variable, Integer> variableMap;
        private List<Rule> rules;

        private List<GroundRule> groundRules;

        public GroundWorker(AtomManager atomManager, GroundRuleStore groundRuleStore,
                Map<Variable, Integer> variableMap, List<Rule> rules) {
            this.atomManager = atomManager;
            this.groundRuleStore = groundRuleStore;
            this.variableMap = variableMap;
            this.rules = rules;

            groundRules = new ArrayList<GroundRule>(rules.size() * PARALLEL_WORKER_QUEUE_SIZE);
        }

        @Override
        public Parallel.Worker<Constant[]> copy() {
            return new GroundWorker(atomManager, groundRuleStore, variableMap, rules);
        }

        @Override
        public void work(int index, Constant[] row) {
            for (Rule rule : rules) {
                GroundRule groundRule = rule.ground(row, variableMap, atomManager);
                if (groundRule != null) {
                    groundRules.add(groundRule);
                }
            }
        }

        /**
         * Dump all the ground rules into the store.
         */
        @Override
        public void batchComplete() {
            groundRuleStore.addAllGroundRules(groundRules);
            groundRules.clear();
        }
    }
}
