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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collection;
import java.util.PriorityQueue;

/**
 * Defines the strategy for exploring the qurey rewrite space.
 */
public abstract class RewriteFringe {
    private static final Logger log = LoggerFactory.getLogger(RewriteFringe.class);

    protected PriorityQueue<RewriteNode> bestNodes;
    protected Collection<RewriteNode> fringe;

    protected RewriteFringe(Collection<RewriteNode> fringe) {
        bestNodes = new PriorityQueue<RewriteNode>();
        this.fringe = fringe;
    }

    public int size() {
        return fringe.size();
    }

    public void clear() {
        bestNodes.clear();
        fringe.clear();
    }

    /**
     * Add a node to the fringe.
     * The node will be check if it is the best one right away because the search through
     * the rewrite space may be time/cost bounded,
     * and we want to check for the best node as soon as we pay the estimation cost (SQL EXPLAIN).
     */
    public void push(RewriteNode node) {
        // Skip invalid rewrites.
        if (node == null) {
            return;
        }

        if (pushInternal(node)) {
            bestNodes.add(node);
            log.trace("   Accepted node.");
        } else {
            log.trace("   Rejected node.");
        }
    }

    public RewriteNode getBestNode() {
        return bestNodes.peek();
    }

    public RewriteNode popBestNode() {
        return bestNodes.poll();
    }

    /**
     * Actually add the node to the fringe.
     * @return true is the node was accepted, false otherwise.
     */
    protected abstract boolean pushInternal(RewriteNode node);

    /**
     * Get the next node for exploration.
     */
    public abstract RewriteNode pop();
}
