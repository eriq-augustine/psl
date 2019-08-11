/**
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

import org.linqs.psl.model.atom.RandomVariableAtom;
import org.linqs.psl.util.RandUtils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Iterate over all the terms from the disk cache.
 * On these non-initial iterations, we will fill the term cache from disk and drain it.
 */
public class DCDStreamingCacheIterator implements DCDStreamingIterator {
    private static final Logger log = LoggerFactory.getLogger(DCDStreamingCacheIterator.class);

    private DCDStreamingTermStore parentStore;
    private Map<Integer, RandomVariableAtom> variables;
    private List<Integer> shuffleMap;

    private List<DCDObjectiveTerm> termCache;
    private List<DCDObjectiveTerm> termPool;

    private ByteBuffer termBuffer;
    private ByteBuffer lagrangeReadBuffer;
    private ByteBuffer lagrangeWriteBuffer;

    private int currentPage;
    private int nextCachedTermIndex;

    private DCDObjectiveTerm nextTerm;

    private boolean shufflePage;

    // When we are reading pages from disk (after the initial round),
    // this list will tell us the order to read them in.
    // This may get shuffled depending on configuration.
    private List<Integer> pageAccessOrder;

    private boolean closed;

    private int numPages;

    private Thread readThread;
    private Thread writeThread;

    public DCDStreamingCacheIterator(
            DCDStreamingTermStore parentStore, Map<Integer, RandomVariableAtom> variables,
            List<DCDObjectiveTerm> termCache, List<DCDObjectiveTerm> termPool,
            ByteBuffer termBuffer, ByteBuffer lagrangeReadBuffer, ByteBuffer lagrangeWriteBuffer,
            boolean shufflePage, List<Integer> shuffleMap, boolean randomizePageAccess,
            int numPages) {
        this.parentStore = parentStore;
        this.variables = variables;
        this.shuffleMap = shuffleMap;

        this.termCache = termCache;
        this.termCache.clear();

        this.termPool = termPool;

        this.termBuffer = termBuffer;
        this.lagrangeReadBuffer = lagrangeReadBuffer;
        this.lagrangeWriteBuffer = lagrangeWriteBuffer;

        nextCachedTermIndex = 0;

        currentPage = -1;
        this.numPages = numPages;

        this.shufflePage = shufflePage;

        pageAccessOrder = new ArrayList<Integer>(numPages);
        for (int i = 0; i < numPages; i++) {
            pageAccessOrder.add(i);
        }

        if (randomizePageAccess) {
            RandUtils.shuffle(pageAccessOrder);
        }

        closed = false;

        // Note that we cannot pre-fetch.
        nextTerm = null;

        readThread = null;
        writeThread = null;

        prefetchPage(currentPage + 1);
    }

    /**
     * Get the next term.
     * It is critical that every call to hasNext be followed by a call to next
     * (as long as hasNext returns true).
     */
    public boolean hasNext() {
        if (nextTerm != null) {
            throw new IllegalStateException("hasNext() was called twice in a row. Call next() directly after hasNext() == true.");
        }

        if (closed) {
            return false;
        }

        nextTerm = fetchNextTerm();
        if (nextTerm == null) {
            close();
            return false;
        }

        return true;
    }

    @Override
    public DCDObjectiveTerm next() {
        if (nextTerm == null) {
            throw new IllegalStateException("Called next() when hasNext() == false (or before the first hasNext() call).");
        }

        DCDObjectiveTerm term = nextTerm;
        nextTerm = null;

        return term;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Get the next term from wherever we need to.
     * We will always settle outstanding pages before trying to get the next term.
     */
    private DCDObjectiveTerm fetchNextTerm() {
        // The cache is exhaused, fill it up.
        if (nextCachedTermIndex >= termCache.size()) {
            // Flush all the lagrange terms.
            flushCache();

            // Check if there is another page, and load it if it exists.
            if (!fetchPage()) {
                // There are no more pages, we are done.
                return null;
            }
        }

        // The page has been verified, so there must be a term waiting.

        DCDObjectiveTerm term = termCache.get(nextCachedTermIndex);
        nextCachedTermIndex++;
        return term;
    }

    /**
     * Fetch the next page.
     * @return true if the next page was fetched and loaded (false when there are no more pages).
     */
    private boolean fetchPage() {
        // Wait for the prefetch to complete.
        waitForPrefetch();

        // Clear the existing page cache.
        termCache.clear();

        currentPage++;
        nextCachedTermIndex = 0;

        if (currentPage >= numPages) {
            // Out of pages.
            return false;
        }

        int numTerms = 0;

        // The buffers are ready from the prefetch,
        // but we still need to convert the data to terms.

        // First get the term size information.
        // We don't actually need to first int (term size) here (only in the prefetch).
        termBuffer.getInt();
        numTerms = termBuffer.getInt();

        // Pull terms out to the pool for reuse.
        for (int i = 0; i < numTerms; i++) {
            DCDObjectiveTerm term = termPool.get(i);
            term.read(termBuffer, lagrangeReadBuffer, variables);
            termCache.add(term);
        }

        if (shufflePage) {
            shuffleMap.clear();
            for (int i = 0; i < termCache.size(); i++) {
                shuffleMap.add(i);
            }

            RandUtils.pairedShuffle(termCache, shuffleMap);
        }

        // Prefetch the next page.
        prefetchPage(currentPage + 1);

        return true;
    }

    private void flushCache() {
        // We don't need to flush if there is nothing to flush.
        if (termCache.size() == 0) {
            return;
        }

        // We will clear the termCache when we fetch a new page, not on flush.
        flushLagrangeCache();
    }

    private void flushLagrangeCache() {
        // Make sure that last write is complete.
        waitForWrite();

        int lagrangeBufferSize = (Float.SIZE / 8) * termCache.size();

        // The buffer has already grown to maximum size in the initial round,
        // no need to reallocate.
        lagrangeWriteBuffer.clear();

        // If this page was picked up from the cache (and not from grounding) and shuffled,
        // then we will need to use the shuffle map to write the lagrange values back in
        // the same order as the terms.
        if (shufflePage) {
            for (int shuffledIndex = 0; shuffledIndex < shuffleMap.size(); shuffledIndex++) {
                int writeIndex = shuffleMap.get(shuffledIndex);
                DCDObjectiveTerm term = termCache.get(shuffledIndex);
                lagrangeWriteBuffer.putFloat(writeIndex * (Float.SIZE / 8), term.getLagrange());
            }
        } else {
            for (DCDObjectiveTerm term : termCache) {
                lagrangeWriteBuffer.putFloat(term.getLagrange());
            }
        }

        asyncWrite(currentPage, lagrangeWriteBuffer, lagrangeBufferSize);
    }

    @Override
    public void close() {
        if (closed) {
            return;
        }
        closed = true;

        flushCache();
        waitForWrite();

        parentStore.cacheIterationComplete();
    }

    private synchronized void prefetchPage(int page) {
        if (readThread != null) {
            throw new IllegalStateException(String.format(
                    "Asked to prefetch a page (%d) when there is already a prefetch in progress.",
                    page));
        }

        if (page >= numPages) {
            return;
        }

        int pageIndex = pageAccessOrder.get(page).intValue();
        String termPagePath = parentStore.getTermPagePath(pageIndex);
        String lagrangePagePath = parentStore.getLagrangePagePath(pageIndex);

        readThread = new ReadWorker(termPagePath, lagrangePagePath);
        readThread.setDaemon(true);
        readThread.start();
    }

    private synchronized void waitForPrefetch() {
        if (readThread == null) {
            return;
        }

        try {
            readThread.join();
        } catch (InterruptedException ex) {
            // Just log and keep going.
            log.warn("Waiting for the prefetch thread was intrrrupted.", ex);
        }

        readThread = null;
    }

    private synchronized void asyncWrite(int page, ByteBuffer buffer, int size) {
        if (writeThread != null) {
            throw new IllegalStateException(String.format(
                    "Asked to write a page (%d) when there is already a write in progress.",
                    page));
        }

        int pageIndex = pageAccessOrder.get(page).intValue();
        String lagrangePagePath = parentStore.getLagrangePagePath(pageIndex);

        writeThread = new WriteWorker(lagrangePagePath, buffer, size);
        writeThread.setDaemon(true);
        writeThread.start();
    }

    private synchronized void waitForWrite() {
        if (writeThread == null) {
            return;
        }

        try {
            writeThread.join();
        } catch (InterruptedException ex) {
            // Just log and keep going.
            log.warn("Waiting for the write thread was intrrrupted.", ex);
        }

        writeThread = null;
    }

    private class ReadWorker extends Thread {
        private String termPagePath;
        private String lagrangePagePath;

        public ReadWorker(String termPagePath, String lagrangePagePath) {
            this.termPagePath = termPagePath;
            this.lagrangePagePath = lagrangePagePath;
        }

        @Override
        public void run() {
            // Prep for the next read.
            // Note that the termBuffer should be at maximum size from the initial round.
            termBuffer.clear();
            lagrangeReadBuffer.clear();

            int termsSize = 0;
            int numTerms = 0;
            int headerSize = (Integer.SIZE / 8) * 2;
            int lagrangesSize = 0;

            try (
                    FileInputStream termStream = new FileInputStream(termPagePath);
                    FileInputStream lagrangeStream = new FileInputStream(lagrangePagePath)) {
                // Read the term size information so we know how much to later read.
                termStream.read(termBuffer.array(), 0, headerSize);

                termsSize = termBuffer.getInt();
                numTerms = termBuffer.getInt();
                lagrangesSize = (Float.SIZE / 8) * numTerms;

                // Now read in all the terms and lagrange values.
                termStream.read(termBuffer.array(), headerSize, termsSize);
                lagrangeStream.read(lagrangeReadBuffer.array(), 0, lagrangesSize);
            } catch (IOException ex) {
                throw new RuntimeException(String.format("Unable to read cache pages: [%s ; %s].", termPagePath, lagrangePagePath), ex);
            }

            // Reset the term buffer so it can read the header again.
            termBuffer.position(0);
        }
    }

    private class WriteWorker extends Thread {
        private String path;
        private ByteBuffer buffer;
        private int size;

        public WriteWorker(String path, ByteBuffer buffer, int size) {
            this.path = path;
            this.buffer = buffer;
            this.size = size;
        }

        @Override
        public void run() {
            try (FileOutputStream stream = new FileOutputStream(path)) {
                stream.write(buffer.array(), 0, size);
            } catch (IOException ex) {
                throw new RuntimeException("Unable to write cache page: " + path, ex);
            }
        }
    }
}
