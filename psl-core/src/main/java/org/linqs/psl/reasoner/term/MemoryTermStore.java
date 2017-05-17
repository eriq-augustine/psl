/*
 * This file is part of the PSL software.
 * Copyright 2011-2015 University of Maryland
 * Copyright 2013-2015 The Regents of the University of California
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
package org.linqs.psl.reasoner.term;

import org.linqs.psl.config.ConfigBundle;
import org.linqs.psl.model.rule.GroundRule;
import org.linqs.psl.model.rule.WeightedGroundRule;

import org.apache.commons.collections4.SetValuedMap;
import org.apache.commons.collections4.multimap.HashSetValuedHashMap;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

public class MemoryTermStore<E extends Term> implements TermStore<E> {
	public static final String CONFIG_PREFIX = "memorytermstore";

	/**
	 * Initial size for the memory store.
	 */
	public static final String INITIAL_SIZE_KEY = CONFIG_PREFIX + ".initialsize";
	public static final int INITIAL_SIZE_DEFAULT = 5000;

	private List<E> store;

	/**
	 * A mapping of ground rule to the term indexes associated with
	 * that ground rule.
	 * This is used for updating weights, so we only track weighted
	 * ground rules and terms.
	 * Note that it could be possible to generate multiple terms
	 * for a single ground rule.
	 */
	private SetValuedMap<WeightedGroundRule, Integer> ruleMapping;

	public MemoryTermStore() {
		this(INITIAL_SIZE_DEFAULT);
	}

	public MemoryTermStore(ConfigBundle config) {
		this(config.getInt(INITIAL_SIZE_KEY, INITIAL_SIZE_DEFAULT));
	}

	public MemoryTermStore(int initialSize) {
		store = new ArrayList<E>(initialSize);
		ruleMapping = new HashSetValuedHashMap<WeightedGroundRule, Integer>(initialSize);
	}

	@Override
	public void add(GroundRule rule, E term) {
		if (rule instanceof WeightedGroundRule && term instanceof WeightedTerm) {
			ruleMapping.put((WeightedGroundRule)rule, new Integer(store.size()));
		}

		store.add(term);
	}

	@Override
	public void clear() {
		store.clear();
		ruleMapping.clear();
	}

	@Override
	public void close() {
		clear();
		store = null;
		ruleMapping = null;
	}

	@Override
	public E get(int index) {
		return store.get(index);
	}

	@Override
	public int size() {
		return store.size();
	}

	@Override
	public Iterator<E> iterator() {
		return store.iterator();
	}

	@Override
	public void updateWeight(WeightedGroundRule rule) {
		for (Integer index : ruleMapping.get(rule)) {
			((WeightedTerm)store.get(index.intValue())).setWeight(rule.getWeight().getWeight());
		}
	}
}