#include "parallel_delta_stepping_v2.h"
#include "delta_stepping.h"
#include <limits>
#include <thread>
#include <mutex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <atomic>


/*
Main differences between the two version :

1. Got rid of all those global helper functions and stuffed everything into a ParallelSolver class.

2. Bucket handling:
   - v1: Used deques and had this findNextBucket thing to skip empty ones
   - v2: Just a vector of vectors, loop through all buckets linearly whihc is simpler

3. Request generation totally changed:
   - v1: Had a separate findRequests() function that built lists
   - v2: Each thread builds its own requests while scanning nodes
   - No more shared request lists = less synchronization headaches

4. The big performance fix:
   - v1: one global mutex for everything terrible idea)
   - v2: Each bucket gets its own mutex
   - Also each thread has its own WorkerData struct to collect stuff locally

5. Relaxation is cleaner:
   - v1: the logic was to find and remove nodes from deques
   - v2: Just collect everything, sort it, pick the best distance per node
   - Don't need to remove nodes from buckets anymore

Basically v2 fixes the fact that global mutex was killing performance. 
Now threads mostly work independently and only lock when updating specific buckets.
*/


const int INF = std::numeric_limits<int>::max();

int findDelta(const Graph& g) {
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

int optimalNumberOfThreads(const Graph& g, int maxThreads) {
    int n = g.size();
    int edge_count = 0;
    for (int u = 0; u < n; u++) {
        edge_count += g.getEdges(u).size();
    }
    if (n < 1000 || edge_count < 10000)
        return 1;
    
    return maxThreads;
}

struct WorkerData {
    std::vector<std::pair<int, int>> light_reqs;
    std::vector<std::pair<int, int>> heavy_reqs;
    WorkerData() = default;
};

class ParallelSolver {
    const Graph& graph;
    int delta;
    int num_workers;
    std::vector<std::atomic<int>> distances;
    std::vector<std::mutex> bucket_locks;
    std::vector<std::vector<int>> buckets;
    
public:
    ParallelSolver(const Graph& g, int d, int threads) : graph(g), delta(d), num_workers(threads), distances(g.size()) {
        // to figure how many buckets we need
        int max_w = 0;
        for (int u = 0; u < graph.size(); u++) {
            for (const Edge& e : graph.getEdges(u)) {
                if (e.weight > max_w) max_w = e.weight;
            }
        }
        
        int bucket_count = (max_w / delta) + 2;
        buckets.resize(bucket_count);
        bucket_locks = std::vector<std::mutex>(bucket_count);
        for (int i = 0; i < graph.size(); i++) {
            distances[i] = INF;
        }
    }
    
    std::vector<int> solve(int start) {
        distances[start] = 0;
        bucket_locks[0].lock();
        buckets[0].push_back(start);
        bucket_locks[0].unlock();
        for (int bucket_idx = 0; bucket_idx < buckets.size(); bucket_idx++) {
            process_bucket(bucket_idx);
        }
        std::vector<int> result(graph.size());
        for (int i = 0; i < graph.size(); i++) {
            result[i] = distances[i].load();
        }
        return result;
    }
    
private:
    void process_bucket(int b_idx) {
        std::vector<int> finished_nodes;
        while (true) {
            std::vector<int> current_batch;
            bucket_locks[b_idx].lock();
            if (buckets[b_idx].empty()) {
                bucket_locks[b_idx].unlock();
                break;
            }
            current_batch = std::move(buckets[b_idx]);
            bucket_locks[b_idx].unlock();
            finished_nodes.insert(finished_nodes.end(), 
                                current_batch.begin(), current_batch.end());
            
            process_light_edges(current_batch, b_idx);
        }
        if (!finished_nodes.empty()) {
            process_heavy_edges(finished_nodes);
        }
    }
    
    void process_light_edges(const std::vector<int>& nodes, int bucket) {
        if (nodes.empty()) return;
        
        std::vector<WorkerData> worker_storage(num_workers);
        std::vector<std::thread> workers;
        
        int per_thread = (nodes.size() + num_workers - 1) / num_workers;
        
        for (int tid = 0; tid < num_workers; tid++) {
            int from = tid * per_thread;
            int to = std::min(from + per_thread, (int)nodes.size());
            
            if (from >= to) continue;
            
            workers.push_back(std::thread([&, tid, from, to]() {
                auto& my_data = worker_storage[tid];
                
                for (int idx = from; idx < to; idx++) {
                    int node = nodes[idx];
                    int d = distances[node].load();
                    for (const Edge& e : graph.getEdges(node)) {
                        if (e.weight <= delta) { 
                            int new_d = d + e.weight;
                            if (new_d < distances[e.dest].load()) {
                                my_data.light_reqs.push_back({e.dest, new_d});
                            }
                        }
                    }
                }
            }));
        }
        for (auto& w : workers) {
            w.join();
        }
        
        do_relaxations(worker_storage, true);
    }
    
    void process_heavy_edges(const std::vector<int>& nodes) {
        if (nodes.empty()) return;
        std::vector<WorkerData> worker_storage(num_workers);
        std::vector<std::thread> workers;
        int per_thread = (nodes.size() + num_workers - 1) / num_workers;
        
        for (int t_id = 0; t_id < num_workers; t_id++) {
            int from = t_id * per_thread;
            int to = std::min(from + per_thread, (int)nodes.size());
            
            if (from >= to) continue;
            
            workers.push_back(std::thread([&, t_id, from, to]() {
                auto& my_data = worker_storage[t_id];
                
                for (int idx = from; idx < to; idx++) {
                    int node = nodes[idx];
                    int d = distances[node].load();
                    
                    for (const Edge& e : graph.getEdges(node)) {
                        if (e.weight > delta) {
                            int new_d = d + e.weight;
                            if (new_d < distances[e.dest].load()) {
                                my_data.heavy_reqs.push_back({e.dest, new_d});
                            }
                        }
                    }
                }
            }));
        }
        
        for (auto& w : workers) {
            w.join();
        }
        
        do_relaxations(worker_storage, false);
    }
    
    void do_relaxations(std::vector<WorkerData>& data, bool is_light) {
        std::vector<std::pair<int, int>> all_reqs;
        
        for (auto& wd : data) {
            auto& reqs = is_light ? wd.light_reqs : wd.heavy_reqs;
            all_reqs.insert(all_reqs.end(), reqs.begin(), reqs.end());
            reqs.clear();
        }
        
        if (all_reqs.empty()) return;
        
        std::sort(all_reqs.begin(), all_reqs.end());
        
        std::vector<std::pair<int, int>> best_reqs;
        int last_node = -1;
        int best_d = INF;
        
        for (const auto& r : all_reqs) {
            if (r.first != last_node) {
                if (last_node >= 0 && best_d < INF) {
                    best_reqs.push_back({last_node, best_d});
                }
                last_node = r.first;
                best_d = r.second;
            } else {
                best_d = std::min(best_d, r.second);
            }
        }
        
        if (last_node >= 0 && best_d < INF) {
            best_reqs.push_back({last_node, best_d});
        }
    
        std::vector<std::thread> workers;
        int per_thread = (best_reqs.size() + num_workers - 1) / num_workers;
        
        for (int tid = 0; tid < num_workers; tid++) {
            int from = tid * per_thread;
            int to = std::min(from + per_thread, (int)best_reqs.size());
            
            if (from >= to) continue;
            
            workers.push_back(std::thread([&, from, to]() {
                for (int i = from; i < to; i++) {
                    int node = best_reqs[i].first;
                    int new_d = best_reqs[i].second;
                    int old_d = distances[node].load();
                    while (new_d < old_d) {
                        if (distances[node].compare_exchange_weak(old_d, new_d)) {
                            int new_b = new_d / delta;
                            if (new_b < buckets.size()) {
                                bucket_locks[new_b].lock();
                                buckets[new_b].push_back(node);
                                bucket_locks[new_b].unlock();
                            }
                            break;
                        }
                    }
                }
            }));
        }
        
        for (auto& w : workers) {
            w.join();
        }
    }
};

std::vector<int> parallelDeltaStepping_v2(const Graph& g, int source, int delta, int numThreads) {
    delta = findDelta(g);
    if (numThreads <= 0) {
        int hw_threads = std::thread::hardware_concurrency();
        numThreads = optimalNumberOfThreads(g, hw_threads);
    }
    if (g.size() < 1000) {
        return deltaStepping(g, source, delta);
    }
    ParallelSolver solver(g, delta, numThreads);
    return solver.solve(source);
}