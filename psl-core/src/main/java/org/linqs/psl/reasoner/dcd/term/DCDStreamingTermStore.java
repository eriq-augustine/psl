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
package org.linqs.psl.reasoner.dcd.term;

import org.linqs.psl.config.Config;
import org.linqs.psl.database.atom.AtomManager;
import org.linqs.psl.model.atom.RandomVariableAtom;
import org.linqs.psl.model.rule.arithmetic.AbstractArithmeticRule;
import org.linqs.psl.model.rule.GroundRule;
import org.linqs.psl.model.rule.Rule;
import org.linqs.psl.model.rule.WeightedRule;
import org.linqs.psl.model.rule.logical.WeightedLogicalRule;
import org.linqs.psl.util.RandUtils;
import org.linqs.psl.util.SystemUtils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.nio.ByteBuffer;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * A term store that iterates over ground queries directly (obviating the GroundRuleStore).
 * Note that the iterators given by this class are meant to be exhaustd (at least the first time).
 * Remember that this class will internally iterate over an unknown number of groundings.
 * So interrupting the iteration can cause the term count to be incorrect.
 */
public class DCDStreamingTermStore implements DCDTermStore {
    private static final Logger log = LoggerFactory.getLogger(DCDStreamingTermStore.class);

    /**
     * Prefix of property keys used by this class.
     */
    public static final String CONFIG_PREFIX = "dcdstreaming";

    /**
     * Where on disk to write term pages.
     */
    public static final String PAGE_LOCATION_KEY = CONFIG_PREFIX + ".pagelocation";
    public static final String PAGE_LOCATION_DEFAULT = SystemUtils.getTempDir("dcd_cache_pages");

    /**
     * The number of terms in a single page.
     */
    public static final String PAGE_SIZE_KEY = CONFIG_PREFIX + ".pagesize";
    public static final int PAGE_SIZE_DEFAULT = 10000;

    /**
     * Whether to shuffle within a page when it is picked up.
     */
    public static final String SHUFFLE_PAGE_KEY = CONFIG_PREFIX + ".shufflepage";
    public static final boolean SHUFFLE_PAGE_DEFAULT = true;

    /**
     * Whether to pick up pages in a random order.
     */
    public static final String RANDOMIZE_PAGE_ACCESS_KEY = CONFIG_PREFIX + ".randomizepageaccess";
    public static final boolean RANDOMIZE_PAGE_ACCESS_DEFAULT = true;

    // How much to over-allocate by.
    public static final double OVERALLOCATION_RATIO = 1.25;

    private List<WeightedLogicalRule> rules;
    private AtomManager atomManager;

    // <Object.hashCode(), RVA>
    private Map<Integer, RandomVariableAtom> variables;

    private boolean initialRound;
    private DCDStreamingIterator activeIterator;
    private int seenTermCount;
    private int numPages;

    private DCDTermGenerator termGenerator;

    private int pageSize;
    private String pageDir;
    private boolean shufflePage;
    private boolean randomizePageAccess;

    /**
     * The IO buffer for terms.
     * This buffer is only written on the first iteration,
     * and contains only components of the terms that do not change.
     */
    private ByteBuffer termBuffer;

    /**
     * The IO buffers for lagrange values.
     * These values change every iteration, and need to be updated.
     */
    private ByteBuffer lagrangeReadBuffer;
    private ByteBuffer lagrangeWriteBuffer;

    /**
     * Terms in the current page.
     * On the initial round, this is filled from DB and flushed to disk.
     * On subsequent rounds, this is filled from disk.
     */
    private List<DCDObjectiveTerm> termCache;

    /**
     * Terms that we will reuse when we start pulling from the cache.
     * This should be a fill page's worth.
     * After the initial round, terms will bounce between here and the term cache.
     */
    private List<DCDObjectiveTerm> termPool;

    /**
     * When we shuffle pages, we need to know how they were shuffled so the lagrange
     * cache can be writtten in the same order.
     * So we will shuffle this list of sequential ints in the same order as the page.
     */
    private List<Integer> shuffleMap;

    public DCDStreamingTermStore(List<Rule> rules, AtomManager atomManager) {
        this.rules = new ArrayList<WeightedLogicalRule>();
        for (Rule rule : rules) {
            if (!rule.isWeighted()) {
                log.warn("DCD does not support hard constraints: " + rule);
                continue;
            }

            // HACK(eriq): This is not actually true,
            //  but I am putting it in place for efficiency reasons.
            if (((WeightedRule)rule).getWeight() < 0.0) {
                log.warn("DCD does not support negative weights: " + rule);
                continue;
            }

            if (rule instanceof AbstractArithmeticRule) {
                log.warn("DCD does not support arithmetic rules: " + rule);
                continue;
            }

            if (!(rule instanceof WeightedLogicalRule)) {
                log.warn("DCD does not support this rule: " + rule);
                continue;
            }

            this.rules.add((WeightedLogicalRule)rule);
        }

        if (rules.size() == 0) {
            throw new IllegalArgumentException("Found no valid rules for DCD.");
        }

        this.atomManager = atomManager;
        termGenerator = new DCDTermGenerator();
        variables = new HashMap<Integer, RandomVariableAtom>();

        initialRound = true;
        activeIterator = null;
        numPages = 0;

        termBuffer = null;
        lagrangeReadBuffer = null;
        lagrangeWriteBuffer = null;

        pageSize = Config.getInt(PAGE_SIZE_KEY, PAGE_SIZE_DEFAULT);
        pageDir = Config.getString(PAGE_LOCATION_KEY, PAGE_LOCATION_DEFAULT);
        shufflePage = Config.getBoolean(SHUFFLE_PAGE_KEY, SHUFFLE_PAGE_DEFAULT);
        randomizePageAccess = Config.getBoolean(RANDOMIZE_PAGE_ACCESS_KEY, RANDOMIZE_PAGE_ACCESS_DEFAULT);
        SystemUtils.recursiveDelete(pageDir);

        if (pageSize <= 1) {
            throw new IllegalArgumentException("Page size is too small.");
        }

        termCache = new ArrayList<DCDObjectiveTerm>(pageSize);
        termPool = new ArrayList<DCDObjectiveTerm>(pageSize);
        shuffleMap = new ArrayList<Integer>(pageSize);

        (new File(pageDir)).mkdirs();
    }

    public boolean isFirstRound() {
        return initialRound;
    }

    @Override
    public int getNumVariables() {
        return variables.size();
    }

    @Override
    public Iterable<RandomVariableAtom> getVariables() {
        return variables.values();
    }

    @Override
    public synchronized RandomVariableAtom createLocalVariable(RandomVariableAtom atom) {
        int key = System.identityHashCode(atom);

        if (variables.containsKey(key)) {
            return atom;
        }

        atom.setValue(RandUtils.nextFloat());
        variables.put(key, atom);

        return atom;
    }

    @Override
    public void ensureVariableCapacity(int capacity) {
        if (capacity == 0) {
            return;
        }

        if (variables.size() == 0) {
            // If there are no variables, then re-allocate the variable storage.
            // The default load factor for Java HashSets is 0.75.
            variables = new HashMap<Integer, RandomVariableAtom>((int)Math.ceil(capacity / 0.75));
        }
    }

    @Override
    public void clear() {
        initialRound = true;
        numPages = 0;

        if (activeIterator != null) {
            activeIterator.close();
            activeIterator = null;
        }

        if (variables != null) {
            variables.clear();
        }

        if (termCache != null) {
            termCache.clear();
        }

        if (termPool != null) {
            termPool.clear();
        }

        SystemUtils.recursiveDelete(pageDir);
    }

    @Override
    public void close() {
        clear();

        if (variables != null) {
            variables = null;
        }

        if (termBuffer != null) {
            termBuffer.clear();
            termBuffer = null;
        }

        if (lagrangeReadBuffer != null) {
            lagrangeReadBuffer.clear();
            lagrangeReadBuffer = null;
        }

        if (lagrangeWriteBuffer != null) {
            lagrangeWriteBuffer.clear();
            lagrangeWriteBuffer = null;
        }

        if (termCache != null) {
            termCache = null;
        }

        if (termPool != null) {
            termPool = null;
        }
    }

    @Override
    public int size() {
        return seenTermCount;
    }

    @Override
    public void add(GroundRule rule, DCDObjectiveTerm term) {
        throw new UnsupportedOperationException();
    }

    @Override
    public DCDObjectiveTerm get(int index) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void ensureCapacity(int capacity) {
        throw new UnsupportedOperationException();
    }

    public String getTermPagePath(int index) {
        return Paths.get(pageDir, String.format("%08d_term.page", index)).toString();
    }

    public String getLagrangePagePath(int index) {
        return Paths.get(pageDir, String.format("%08d_lagrange.page", index)).toString();
    }

    /**
     * A callback for the initial round iterator.
     * The ByterBuffers are here because of possible reallocation.
     */
    public void initialIterationComplete(int termCount, int numPages, ByteBuffer termBuffer, ByteBuffer lagrangeWriteBuffer) {
        seenTermCount = termCount;
        this.numPages = numPages;
        this.termBuffer = termBuffer;
        this.lagrangeWriteBuffer = lagrangeWriteBuffer;

        // Now allocate the lagrange read buffer bassed on the size of the write buffer.
        lagrangeReadBuffer = ByteBuffer.allocate(lagrangeWriteBuffer.capacity());

        initialRound = false;
        activeIterator = null;
    }

    /**
     * A callback for the non-initial round iterator.
     */
    public void cacheIterationComplete() {
        activeIterator = null;
    }

    @Override
    public Iterator<DCDObjectiveTerm> iterator() {
        if (activeIterator != null) {
            throw new IllegalStateException("Iterator already exists for this DCDTermStore. Exhaust the iterator first.");
        }

        if (initialRound) {
            activeIterator = new DCDStreamingInitialRoundIterator(
                    this, rules, atomManager, termGenerator,
                    termCache, termPool, termBuffer, lagrangeWriteBuffer, pageSize);
        } else {
            activeIterator = new DCDStreamingCacheIterator(
                    this, variables, termCache, termPool,
                    termBuffer, lagrangeReadBuffer, lagrangeWriteBuffer, shufflePage, shuffleMap, randomizePageAccess, numPages);
        }

        return activeIterator;
    }
}
