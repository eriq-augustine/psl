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
package org.linqs.psl.grounding.rewrite;

import org.linqs.psl.model.formula.Formula;
import org.linqs.psl.util.MathUtils;

/**
 * A container for rewrite search nodes.
 */
public class RewriteNode implements Comparable<RewriteNode> {
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

    @Override
    public String toString() {
        return String.format("{Atom Bits: %d, Optimistic: %f, Pessimistic: %f, Formula: %s}",
                atomsBitSet, optimisticCost, pessimisticCost, formula);
    }

    @Override
    public int compareTo(RewriteNode other) {
        if (other == null) {
            return 1;
        }

        return MathUtils.compare(this.optimisticCost, other.optimisticCost);
    }

    @Override
    public boolean equals(Object other) {
        if (other == null || !(other instanceof RewriteNode)) {
            return false;
        }

        if (other == this) {
            return true;
        }

        RewriteNode otherNode = (RewriteNode)other;
        return (otherNode.atomsBitSet == this.atomsBitSet) && (otherNode.formula.equals(this.formula));
    }
}
