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
package org.linqs.psl.util;

import org.linqs.psl.config.Config;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;

/**
 * Utilities to run operations in parallel.
 * The threads will be started up on the first call (or initPool()), and not shut down until the JVM shuts down.
 * Since the thread pool (and CPU) is shared, only one task may be run in parallel at a time.
 */
public final class Parallel {
    private static final Logger log = LoggerFactory.getLogger(Parallel.class);

    public static final String CONFIG_PREFIX = "parallel";

    public static final String NUM_THREADS_KEY = CONFIG_PREFIX + ".numthreads";
    public static final int NUM_THREADS_DEFAULT = Runtime.getRuntime().availableProcessors();

    private static boolean initialized = false;

    // Defer assignment until a request is actually made to let the config get initialized.
    // Don't refer to directly, use getNumThreads() instead.
    private static int numThreads = -1;

    // The number of jobs a single worker can have, managed internally.
    private static final int DEFAULT_WORKER_QUEUE_SIZE = 1;

    // Block putting work intot he pool until there are workers ready.
    private static BlockingQueue<Worker<?>> workerQueue;

    // Keep all the workers somewhere we can reference them.
    private static List<Worker<?>> allWorkers;

    private static ExecutorService pool;

    // Timing numbers used in all parallel methods.
    // Allocating here instead of per-invocation allowes easier method decomposition.
    // (Since we don't have to worry about allocating a wrapper to pass them around.)
    private static long iterations = 0;
    private static long parentWaitTimeMS = 0;
    private static long workerWaitTimeMS = 0;
    private static long workerWorkTimeMS = 0;

    /**
     * Objects that are specific to each thread.
     */
    private static Map<Thread, Map<String, Object>> threadObjects = new ConcurrentHashMap<Thread, Map<String, Object>>();

    // Static only.
    private Parallel() {}

    public synchronized static int getNumThreads() {
        if (numThreads == -1) {
            numThreads = Config.getInt(NUM_THREADS_KEY, NUM_THREADS_DEFAULT);
        }

        return numThreads;
    }

    public static boolean hasThreadObject(String key) {
        if (!threadObjects.containsKey(Thread.currentThread())) {
            threadObjects.put(Thread.currentThread(), new HashMap<String, Object>());
        }

        return threadObjects.get(Thread.currentThread()).containsKey(key);
    }

    public static Object getThreadObject(String key) {
        if (!threadObjects.containsKey(Thread.currentThread())) {
            threadObjects.put(Thread.currentThread(), new HashMap<String, Object>());
        }

        return threadObjects.get(Thread.currentThread()).get(key);
    }

    public static void putThreadObject(String key, Object value) {
        if (!threadObjects.containsKey(Thread.currentThread())) {
            threadObjects.put(Thread.currentThread(), new HashMap<String, Object>());
        }

        threadObjects.get(Thread.currentThread()).put(key, value);
    }

    /**
     * Count and call a worker with each number in [start, end).
     * Inclusive with start, exclusive with end.
     * The caller is trusted to provide appropriate numbers.
     */
    public synchronized static RunTimings count(int start, int end, int increment, int workerQueueSize, Worker<Integer> baseWorker) {
        initWorkers(baseWorker, workerQueueSize);
        RunTimings timings = countInternal(start, end, increment, workerQueueSize);
        cleanupWorkers();

        return timings;
    }

    /**
     * Convenience count() that increments by 1.
     */
    public static RunTimings count(int start, int end, Worker<Integer> baseWorker) {
        return count(start, end, 1, DEFAULT_WORKER_QUEUE_SIZE, baseWorker);
    }

    /**
     * Convenience count() that increments by 1.
     */
    public static RunTimings count(int start, int end, int workerQueueSize, Worker<Integer> baseWorker) {
        return count(start, end, 1, workerQueueSize, baseWorker);
    }

    /**
     * Convenience count() that starts at 0 and increments by 1.
     */
    public static RunTimings count(int end, Worker<Integer> baseWorker) {
        return count(0, end, 1, DEFAULT_WORKER_QUEUE_SIZE, baseWorker);
    }

    private static RunTimings countInternal(int start, int end, int increment, int workerQueueSize) {
        resetTiming();

        int number = start;

        while (number < end) {
            Worker<?> worker = fetchOpenWorker();

            @SuppressWarnings("unchecked")
            Worker<Integer> intWorker = (Worker<Integer>)worker;

            for (int i = 0; i < workerQueueSize; i++) {
                if (number >= end) {
                    break;
                }

                intWorker.addWork(number, new Integer(number));
                number += increment;
            }

            if (intWorker.workQueuedCount() > 0) {
                pool.execute(intWorker);
            }
        }

        flushWorkerQueue();

        return new RunTimings();
    }

    /**
     * Invoke a worker once for each item.
     */
    public synchronized static <T> RunTimings foreach(Iterable<T> work, int  workerQueueSize, Worker<T> baseWorker) {
        initWorkers(baseWorker, workerQueueSize);
        RunTimings timings = foreachInternal(work, workerQueueSize);
        cleanupWorkers();

        return timings;
    }

    public static <T> RunTimings foreach(Iterable<T> work, Worker<T> baseWorker) {
        int workerQueueSize = DEFAULT_WORKER_QUEUE_SIZE;
        if (work instanceof List) {
            workerQueueSize = Math.max(1, ((List)work).size() / getNumThreads());
        }

        return foreach(work, workerQueueSize, baseWorker);
    }

    public static <T> RunTimings foreach(Iterator<T> work, Worker<T> baseWorker) {
        return foreach(IteratorUtils.newIterable(work), baseWorker);
    }

    public static <T> RunTimings foreach(Iterator<T> work, int workerQueueSize, Worker<T> baseWorker) {
        return foreach(IteratorUtils.newIterable(work), workerQueueSize, baseWorker);
    }

    private static <T> RunTimings foreachInternal(Iterable<T> work, int workerQueueSize) {
        resetTiming();

        int count = 0;
        Iterator<T> workIterator = work.iterator();

        while (workIterator.hasNext()) {
            Worker<?> worker = fetchOpenWorker();

            @SuppressWarnings("unchecked")
            Worker<T> typedWorker = (Worker<T>)worker;

            for (int i = 0; i < workerQueueSize; i++) {
                if (!workIterator.hasNext()) {
                    break;
                }

                typedWorker.addWork(count, workIterator.next());
                count++;
            }

            if (typedWorker.workQueuedCount() > 0) {
                pool.execute(typedWorker);
            }
        }

        flushWorkerQueue();

        return new RunTimings();
    }

    /**
     * Init the thread pool and supporting structures.
     */
    public static synchronized void initPool() {
        if (initialized) {
            return;
        }

        // We can use an unbounded queue (no initial size given) since the parent
        // thread is disciplined when giving out work.
        workerQueue = new LinkedBlockingQueue<Worker<?>>();
        allWorkers = new ArrayList<Worker<?>>(getNumThreads());

        // We will make all the threads daemons, so the JVM shutdown will not be held up.
        pool = Executors.newFixedThreadPool(getNumThreads(), new DaemonThreadFactory());

        // Close the pool only at JVM shutdown.
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                Parallel.shutdown();
            }
        });

        initialized = true;
    }

    /**
     * Reset internal timings.
     */
    private static void resetTiming() {
        iterations = 0;
        parentWaitTimeMS = 0;
        workerWaitTimeMS = 0;
        workerWorkTimeMS = 0;
    }

    private static Worker<?> fetchOpenWorker() {
        Worker<?> worker = null;
        try {
            // Will block if no workers are ready.
            long time = System.currentTimeMillis();
            worker = workerQueue.take();
            parentWaitTimeMS += (System.currentTimeMillis() - time);
            iterations++;
        } catch (InterruptedException ex) {
            throw new RuntimeException("Interrupted waiting for worker.");
        }

        if (worker.getException() != null) {
            throw new RuntimeException("Exception on worker.", worker.getException());
        }

        return worker;
    }

    /**
     * As workers finish, they will be added to the queue.
     * We can wait for all the workers by emptying out the queue.
     */
    private static void flushWorkerQueue() {
        for (int i = 0; i < getNumThreads(); i++) {
            try {
                long time = System.currentTimeMillis();
                Worker<?> worker = workerQueue.take();
                parentWaitTimeMS += (System.currentTimeMillis() - time);

                workerWaitTimeMS += worker.getWaitTime();
                workerWorkTimeMS += worker.getWorkTime();

                if (worker.getException() != null) {
                    throw new RuntimeException("Exception on worker.", worker.getException());
                }
            } catch (InterruptedException ex) {
                throw new RuntimeException("Interrupted waiting for worker (" + i + ").");
            }
        }
    }

    /**
     * Always the first thing called when setting up to run a task in parallel.
     */
    private static <T> void initWorkers(Worker<T> baseWorker, int workerQueueSize) {
        initPool();

        if (workerQueueSize <= 0) {
            throw new RuntimeException("Worker queue size must be positive, got: " + workerQueueSize + ".");
        }

        workerQueue.clear();
        allWorkers.clear();

        for (int i = 0; i < getNumThreads(); i++) {
            Worker<T> worker = null;

            // The base worker goes in last so we won't call copy() after init().
            if (i == getNumThreads() - 1) {
                worker = baseWorker;
            } else {
                worker = baseWorker.copy();
            }

            worker.init(i, workerQueueSize);

            allWorkers.add(worker);
            workerQueue.add(worker);
        }
    }

    private static void cleanupWorkers() {
        for (Worker<?> worker : allWorkers) {
            worker.close();
        }

        allWorkers.clear();
        workerQueue.clear();
    }

     private static void shutdown() {
        cleanupWorkers();

        try {
            pool.shutdownNow();
            pool.awaitTermination(10, TimeUnit.SECONDS);
        } catch (InterruptedException ex) {
            // Do nothing, we are shutting down anyways.
        }

        workerQueue = null;
        allWorkers = null;
        pool = null;
     }

    /**
     * Signal that a worker is done and ready for more work.
     */
    private static void freeWorker(Worker<?> worker) {
        workerQueue.add(worker);
    }

    /**
     * Extend this class for any work.
     * Default implmentation are provided for all non-abstract, non-final methods.
     */
    public static abstract class Worker<T> implements Runnable, Cloneable {
        protected int id;

        private long waitTimeMS;
        private long workTimeMS;
        /**
         * The queued items to work on.
         * We will be very careful in this class to avoid unnecessary allocations while
         * still maintining a distinct queue for each worker.
         */
        private List<Integer> indexes;
        private List<T> items;

        private Exception exception;

        public Worker() {
            this.id = -1;
            this.waitTimeMS = 0;
            this.workTimeMS = 0;
            this.indexes = new ArrayList<Integer>(DEFAULT_WORKER_QUEUE_SIZE);
            this.items = new ArrayList<T>(DEFAULT_WORKER_QUEUE_SIZE);
            this.exception = null;
        }

        /**
         * Cleanup anything.
         * Called after all work has been complete and it is time to clean up.
         */
        public void close() {}

        /**
         * Cleanup anything between batches.
         * Called after a batch of work has been complete.
         * The size of the bath is determined by the set worker queue size.
         * but may be less.
         */
        public void batchComplete() {}

        @Override
        protected final Object clone() {
            Worker<T> newWorker;

            try {
                @SuppressWarnings("unchecked")
                Worker<T> tempWorker = (Worker<T>)super.clone();
                newWorker = tempWorker;
            } catch (CloneNotSupportedException ex) {
                throw new RuntimeException("Either implement copy(), or support clone() for Workers.", ex);
            }

            // Explicitly deep copy the work queue.
            newWorker.indexes = new ArrayList<Integer>(DEFAULT_WORKER_QUEUE_SIZE);
            newWorker.items = new ArrayList<T>(DEFAULT_WORKER_QUEUE_SIZE);

            return newWorker;
        }

        /**
         * Make a deep copy of this worker.
         * Called when the manager is getting the correct number of workers ready.
         * Children can override, just make sure that either clone() is used or the constructor is called.
         */
        public Worker<T> copy() {
            @SuppressWarnings("unchecked")
            Worker<T> worker = (Worker<T>)this.clone();
            return worker;
        }

        /**
         * Called before any work is given.
         * The id will be unique to this worker for this batch of work.
         */
        public final void init(int id, int capacity) {
            this.id = id;

            ((ArrayList<Integer>)indexes).ensureCapacity(capacity);
            ((ArrayList<T>)items).ensureCapacity(capacity);
        }

        public final void clearException() {
            exception = null;
        }

        public final Exception getException() {
            return exception;
        }

        public final long getWaitTime() {
            return waitTimeMS;
        }

        public final long getWorkTime() {
            return workTimeMS;
        }

        public final int workQueuedCount() {
            return items.size();
        }

        @Override
        public String toString() {
            return getClass().getName() + "::" + id;
        }

        @Override
        public final void run() {
            try {
                if (indexes.size() == 0) {
                    log.warn("Called run() without first calling addWork().");
                    return;
                }

                for (int i = 0; i < items.size(); i++) {
                    long time = System.currentTimeMillis();
                    work(indexes.get(i), items.get(i));
                    workTimeMS += (System.currentTimeMillis() - time);
                }

                batchComplete();
            } catch (Exception ex) {
                log.warn("Caught exception on worker: {}", id);
                exception = ex;
            } finally {
                indexes.clear();
                items.clear();

                long time = System.currentTimeMillis();
                freeWorker(this);
                waitTimeMS += (System.currentTimeMillis() - time);
            }
        }

        public final void addWork(int index, T item) {
            this.indexes.add(index);
            this.items.add(item);
        }

        /**
         * Do the actual work.
         * The index is the item's index in the collection.
         */
        public abstract void work(int index, T item);
    }

    private static class DaemonThreadFactory implements ThreadFactory {
        private ThreadFactory defaultThreadFactory;

        public DaemonThreadFactory() {
            this.defaultThreadFactory = Executors.defaultThreadFactory();
        }

        @Override
        public Thread newThread(Runnable r) {
            Thread thread = defaultThreadFactory.newThread(r);
            thread.setDaemon(true);
            return thread;
        }
    }

    public static class RunTimings {
        public final long iterations;
        public final long parentWaitTimeMS;
        public final long workerWaitTimeMS;
        public final long workerWorkTimeMS;

        public RunTimings() {
            this.iterations = Parallel.iterations;
            this.parentWaitTimeMS = Parallel.parentWaitTimeMS;
            this.workerWaitTimeMS = Parallel.workerWaitTimeMS;
            this.workerWorkTimeMS = Parallel.workerWorkTimeMS;
        }

        public String toString() {
            return String.format("Iterations: %d, Parent Wait Time: %d, Worker Wait Time: %d, Worker Work Time: %d",
                    this.iterations, this.parentWaitTimeMS, this.workerWaitTimeMS, this.workerWorkTimeMS);
        }
    }
}
