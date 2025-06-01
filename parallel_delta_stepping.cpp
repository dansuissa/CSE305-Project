#include "parallel_delta_stepping.h"
#include <limits>
#include <thread>
#include <mutex>
#include <vector>
#include <cmath>
#include <deque>
#include <algorithm>
#include <atomic>

const int INF = std::numeric_limits<int>::max();

using EdgePair = std::pair<int, int>;

/* Finding the optimal delta  was tricky since we didn't follow Meyer's approach exactly. 
   Used LLM to help figure out a good value for delta selection. */

/* All the core parallelization logic and implementaton is our own work: thread management, synchronization strategy and so on */

static int findDelta(const Graph& g) {
    int n = g.size();
    int min_w = INF, max_w = 0;
    for (int u = 0; u < n; u++) {
        for (const Edge& e : g.getEdges(u)) {
            if (e.weight > 0 && e.weight < min_w) 
                min_w = e.weight;
            if (e.weight > max_w) 
                max_w = e.weight;
        }
    }
    
    if (min_w == INF) return 1;
    if (n < 100) return min_w;
    int total_edges = 0;
    for (int u = 0; u < n; u++) {
        total_edges += g.getEdges(u).size();
    }
    int avg_degree = (total_edges + n - 1) / n;
    int delta = std::max(1, min_w * std::min(10, avg_degree));
    return std::min(delta, max_w / 10);
}

static int optimalNumberOfThreads(const Graph& g, int hw_threads) {
    int m = 0;
    for (int u = 0; u < g.size(); ++u) m += g.getEdges(u).size();
    int n = g.size();
    if (n < 1000 || m < 5000) return std::min(2, hw_threads);
    if (n < 10000 || m < 50000) return std::min(std::max(2, hw_threads / 2), hw_threads);
    return hw_threads;
}

static int calculateNumBuckets(const Graph& g, int delta) {
    int max_w = 0;
    for (int u = 0; u < g.size(); ++u)
        for (const Edge& e : g.getEdges(u))
            if (e.weight > max_w) max_w = e.weight;
    return max_w;
}

static int findNextBucket(const std::vector<std::deque<int>>& buckets, int start, int maxBuckets) {
    for (int i = start; i < maxBuckets; ++i) if (!buckets[i].empty()) return i;
    return maxBuckets;
}

static void do_relaxations(const std::vector<EdgePair>& req, int L, int R, std::mutex& mtx, std::vector<std::atomic<int>>& dist, std::vector<std::deque<int>>& buckets, int delta) {
    std::vector<EdgePair> pending; pending.reserve(R - L);
    for (int i = L; i < R && i < (int)req.size(); ++i) {
        int v = req[i].first; int new_d = req[i].second;
        if (new_d < dist[v].load()) pending.emplace_back(v, new_d);
    }
    if (pending.empty()) return;
    std::lock_guard<std::mutex> g(mtx);
    for (auto& p : pending) {
        int v = p.first; int new_d = p.second;
        if (new_d < dist[v].load()) {
            int old_d = dist[v].load();
            if (old_d != INF) {
                int old_bucket = old_d / delta;
                if (old_bucket < (int)buckets.size()) {
                    auto& B = buckets[old_bucket];
                    auto it = std::find(B.begin(), B.end(), v);
                    if (it != B.end()) { if (it != B.end() - 1) std::swap(*it, B.back()); B.pop_back(); }
                }
            }
            int new_bucket = new_d / delta;
            if (new_bucket >= (int)buckets.size()) buckets.resize(new_bucket + 1);
            buckets[new_bucket].push_back(v);
            dist[v].store(new_d);
        }
    }
}

static void findRequests(std::vector<EdgePair>& out, const std::vector<int>& nodes, bool lightEdges, const Graph& g, const std::vector<std::atomic<int>>& dist, int delta) {
    out.clear();
    for (int v : nodes) {
        for (const Edge& e : g.getEdges(v)) { 
            if ((lightEdges && e.weight <= delta) || (!lightEdges && e.weight > delta)) {
                int w = e.dest; int new_d = dist[v].load() + e.weight;
                if (new_d < dist[w].load()) out.emplace_back(w, new_d);
            }
        }
    }
}

std::vector<int> parallelDeltaStepping(const Graph& g, int src, int delta, int numThreads) {
    if (delta <= 0) delta = findDelta(g);
    int hw_threads = std::thread::hardware_concurrency(); if (hw_threads == 0) hw_threads = 4;
    int threads = (numThreads > 0) ? std::min(numThreads, hw_threads) : optimalNumberOfThreads(g, hw_threads);
    if (g.size() < 500) threads = 1;

    int bucketCount = calculateNumBuckets(g, delta);
    std::vector<std::deque<int>> buckets(bucketCount);
    std::vector<std::atomic<int>> dist(g.size()); for (auto& d : dist) d.store(INF);
    std::mutex mtx;
    dist[src].store(0); buckets[0].push_back(src);

    auto parallelRelax = [&](std::vector<EdgePair>& req) {
        if (req.empty()) return;
        int chunk = (req.size() + threads - 1) / threads;
        std::vector<std::thread> pool;
        for (int t = 0; t < threads - 1; ++t) {
            int L = t * chunk; int R = std::min<int>((t + 1) * chunk, req.size());
            if (L < R) pool.emplace_back(do_relaxations, std::cref(req), L, R, std::ref(mtx), std::ref(dist), std::ref(buckets), delta);
        }
        int L0 = (threads - 1) * chunk; if (L0 < (int)req.size()) do_relaxations(req, L0, req.size(), mtx, dist, buckets, delta);
        for (auto& th : pool) th.join();
    };

        int k = findNextBucket(buckets, 0, buckets.size());
    while (k < (int)buckets.size()) {
        std::vector<int> settled, R;
        {
            std::lock_guard<std::mutex> lk(mtx);
            R.assign(buckets[k].begin(), buckets[k].end());
            buckets[k].clear();
        }
        while (!R.empty()) {
            settled.insert(settled.end(), R.begin(), R.end());
            std::vector<EdgePair> light;
            findRequests(light, R, true, g, dist, delta);
            parallelRelax(light);
            std::lock_guard<std::mutex> lk(mtx);
            R.assign(buckets[k].begin(), buckets[k].end());
            buckets[k].clear();
        }
        std::vector<EdgePair> heavy;
        findRequests(heavy, settled, false, g, dist, delta);
        parallelRelax(heavy);
        k = findNextBucket(buckets, k + 1, buckets.size());
    }

    std::vector<int> res(g.size());
    for (int i = 0; i < g.size(); ++i) res[i] = dist[i].load();
    return res;
}