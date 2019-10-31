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

import java.util.PriorityQueue;

public class UCSRewriteFringe extends RewriteFringe {
    public UCSRewriteFringe() {
        // Use a priority queue (with PriorityQueue's default size of 11.
        super(new PriorityQueue<RewriteNode>(11));
    }

    protected boolean pushInternal(RewriteNode node) {
        ((PriorityQueue<RewriteNode>)fringe).add(node);
        return true;
    }

    public RewriteNode pop() {
        return ((PriorityQueue<RewriteNode>)fringe).poll();
    }
}
