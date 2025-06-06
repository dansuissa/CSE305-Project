#include "graph.h"
#include "delta_stepping.h"
#include "dijkstra.h"
#include "parallel_delta_stepping.h"
#include "parallel_delta_stepping_v2.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <thread>
#include <fstream>
#include <cmath>
#include <random>

const int DEFAULT_DELTA = 50;
const int MEASUREMENT_TRIALS = 5;

double measureMs(std::function<void()> func) {
    auto start = std::chrono::steady_clock::now();
    func();
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double getAverageTime(std::function<void()> func, int trials = MEASUREMENT_TRIALS) {
    double totalTime = 0.0;
    for (int i = 0; i < trials; i++) {
        totalTime += measureMs(func);
    }
    return totalTime / trials;
}

void writeResultsToFile(const std::string& filename, const std::vector<std::vector<std::string>>& data) {
    std::ofstream file(filename);
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); i++) {
            file << row[i];
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}

bool runCorrectnessTests() {
    std::cout << "Running correctness tests...\n";
    int threadCount = std::thread::hardware_concurrency();
//    std::cout << "[DEBUG] runCorrectnessTests: hardware_concurrency() returned " 
//               << threadCount << " threads\n";
 
    std::vector<std::tuple<std::string, int, int>> testConfigs = {
        {"random_sparse", 100, 2},
        {"random_dense", 500, 10},
        {"grid_small", 900, 0},
        {"grid_large", 10000, 0},
        {"scale_free", 1000, 3},
        {"scale_free_large", 5000, 5},
        {"small_world", 1000, 6},
        {"random_large", 10000, 4}
    };
    for (const auto& [type, size, edge_factor] : testConfigs) {
        Graph testGraph(1);
        if (type.find("random") != std::string::npos) {
            testGraph = Graph::randomGraph(size, size * edge_factor, 100);
        } else if (type.find("grid") != std::string::npos) {
            int side = (int)sqrt(size);
            testGraph = Graph::gridGraph(side, side, 50);
        } else if (type.find("scale_free") != std::string::npos) {
            testGraph = Graph::scaleFreeGraph(size, edge_factor, 100);
        } else if (type.find("small_world") != std::string::npos) {
            testGraph = Graph::smallWorldGraph(size, edge_factor, 0.3, 50);
        }
        auto dijkstraResult = dijkstra(testGraph, 0);
        auto deltaResult = deltaStepping(testGraph, 0, DEFAULT_DELTA);
        auto parallelV1Result = parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount);
        auto parallelV2Result = parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount);
        bool allPass = true;
        if (dijkstraResult != deltaResult) {
            allPass = false;
            std::cout << "Test failed: " << type << " (size=" << size << ") - delta-stepping\n";
            for (int i = 0; i < dijkstraResult.size(); i++) {
                if (dijkstraResult[i] != deltaResult[i]) {
                    std::cout << "  First mismatch at node " << i
                              << ": dijkstra=" << dijkstraResult[i]
                              << ", delta=" << deltaResult[i] << "\n";
                    break;
                }
            }
        }
        // Check parallel v1
        if (dijkstraResult != parallelV1Result) {
            allPass = false;
            std::cout << "Test failed: " << type << " (size=" << size << ") - parallel_v1\n";
            for (int i = 0; i < dijkstraResult.size(); i++) {
                if (dijkstraResult[i] != parallelV1Result[i]) {
                    std::cout << "  First mismatch at node " << i
                              << ": dijkstra=" << dijkstraResult[i]
                              << ", parallel_v1=" << parallelV1Result[i] << "\n";
                    break;
                }
            }
        }
        // Check parallel v2
        if (dijkstraResult != parallelV2Result) {
            allPass = false;
            std::cout << "Test failed: " << type << " (size=" << size << ") - parallel_v2\n";
            for (int i = 0; i < dijkstraResult.size(); i++) {
                if (dijkstraResult[i] != parallelV2Result[i]) {
                    std::cout << "  First mismatch at node " << i
                              << ": dijkstra=" << dijkstraResult[i]
                              << ", parallel_v2=" << parallelV2Result[i] << "\n";
                    break;
                }
            }
        }
        
        if (!allPass) {
            return false;
        }
        std::cout << "OK " << type << " (n=" << size << ") - all algorithms match\n";
    }
    std::cout << "\nAll correctness tests passed! \n\n";
    return true;
}

void analyzeGraphTypes() {
    std::cout << "Analyzing different graph types...\n";
    
    std::vector<std::vector<std::string>> results;
    results.push_back({"graph_type", "size", "edges", "dijkstra_time", "delta_time", "parallel_v1_time", "parallel_v2_time", "seq_speedup", "par_v1_speedup", "par_v2_speedup"});
    
    std::vector<int> testSizes = {5000, 15000,30000, 50000};
    std::vector<std::string> graphTypes = {"random", "grid", "scale_free", "small_world"};
    
    int threadCount = std::thread::hardware_concurrency();
    std::cout << threadCount;
    // if (threadCount == 0) threadCount = 1;
    
    for (const std::string& type : graphTypes) {
        std::cout << "Testing " << type << " graphs...\n";
        
        for (int n : testSizes) {
            Graph testGraph(1);
            int actualEdges = 0;
            
            if (type == "random") {
                int edgeCount = n * 3;
                testGraph = Graph::randomGraph(n, edgeCount, 100);
                actualEdges = edgeCount;
            } else if (type == "grid") {
                int side = (int)sqrt(n);
                testGraph = Graph::gridGraph(side, side, 50);
                actualEdges = 2 * (side - 1) * side;
            } else if (type == "scale_free") {
                testGraph = Graph::scaleFreeGraph(n, 3, 100);
                for (int i = 0; i < n; i++) {
                    actualEdges += testGraph.getEdges(i).size();
                }
            } else if (type == "small_world") {
                testGraph = Graph::smallWorldGraph(n, 6, 0.3, 50);
                for (int i = 0; i < n; i++) {
                    actualEdges += testGraph.getEdges(i).size();
                }
            }
            
            if (testGraph.size() < 50) continue;
            
            double dijkstraTime = getAverageTime([&]() { dijkstra(testGraph, 0); });
            double deltaTime = getAverageTime([&]() { deltaStepping(testGraph, 0, DEFAULT_DELTA); });
            double parallelV1Time = getAverageTime([&]() { parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount); });
            double parallelV2Time = getAverageTime([&]() { parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount); });
            
            double sequentialSpeedup = dijkstraTime / deltaTime;
            double parallelV1Speedup = deltaTime / parallelV1Time;
            double parallelV2Speedup = deltaTime / parallelV2Time;
            
            results.push_back({
                type,
                std::to_string(testGraph.size()),
                std::to_string(actualEdges),
                std::to_string(dijkstraTime),
                std::to_string(deltaTime),
                std::to_string(parallelV1Time),
                std::to_string(parallelV2Time),
                std::to_string(sequentialSpeedup),
                std::to_string(parallelV1Speedup),
                std::to_string(parallelV2Speedup)
            });
            
            std::cout << "  Size " << testGraph.size() << ": v1=" << std::fixed << std::setprecision(2) << parallelV1Speedup 
                      << "x, v2=" << parallelV2Speedup << "x speedup\n";
        }
    }
    
    writeResultsToFile("graph_type_analysis.csv", results);
    std::cout << "Graph type results saved\n\n";
}

void analyzeDensityEffects() {
    std::cout << "Analyzing density effects...\n";
    
    std::vector<std::vector<std::string>> results;
    results.push_back({"size", "density", "edges", "dijkstra_time", "delta_time", "parallel_v1_time", "parallel_v2_time", "seq_speedup", "par_v1_speedup", "par_v2_speedup"});
    
    std::vector<int> testSizes = {2000, 5000, 10000 };
    std::vector<double> densityLevels = {0.005, 0.02, 0.1, 0.2};
    
    int threadCount = std::thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 1;
    
    for (int n : testSizes) {
        std::cout << "Testing size " << n << "...\n";
        
        for (double density : densityLevels) {
            int edgeCount = (int)(density * n * (n - 1) / 2);
            if (edgeCount < n - 1) edgeCount = n - 1;
            
            Graph testGraph = Graph::randomGraph(n, edgeCount, 100);
            
            double dijkstraTime = getAverageTime([&]() { dijkstra(testGraph, 0); });
            double deltaTime = getAverageTime([&]() { deltaStepping(testGraph, 0, DEFAULT_DELTA); });
            double parallelV1Time = getAverageTime([&]() { parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount); });
            double parallelV2Time = getAverageTime([&]() { parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount); });
            
            double sequentialSpeedup = dijkstraTime / deltaTime;
            double parallelV1Speedup = deltaTime / parallelV1Time;
            double parallelV2Speedup = deltaTime / parallelV2Time;
            
            results.push_back({
                std::to_string(n),
                std::to_string(density),
                std::to_string(edgeCount),
                std::to_string(dijkstraTime),
                std::to_string(deltaTime),
                std::to_string(parallelV1Time),
                std::to_string(parallelV2Time),
                std::to_string(sequentialSpeedup),
                std::to_string(parallelV1Speedup),
                std::to_string(parallelV2Speedup)
            });
            
            std::cout << "  Density " << std::fixed << std::setprecision(3) << density 
                      << ": v1=" << std::setprecision(2) << parallelV1Speedup 
                      << "x, v2=" << parallelV2Speedup << "x speedup\n";
        }
    }
    
    writeResultsToFile("density_analysis.csv", results);
    std::cout << "Density analysis results saved\n\n";
}

void analyzeWeightDistributions() {
    std::cout << "Analyzing weight distributions...\n";
    
    std::vector<std::vector<std::string>> results;
    results.push_back({"distribution", "size", "edges", "dijkstra_time", "delta_time", "parallel_v1_time", "parallel_v2_time", "seq_speedup", "par_v1_speedup", "par_v2_speedup"});
    
    std::vector<int> testSizes = {5000, 10000, 30000, 50000};
    std::vector<std::string> distributions = {"uniform", "exponential", "normal", "bimodal"};
    
    int threadCount = std::thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 1;
    
    for (const std::string& dist : distributions) {
        std::cout << "Testing " << dist << " distribution...\n";
        
        for (int n : testSizes) {
            int edgeCount = n * 3;
            Graph testGraph = Graph::randomGraphWithWeights(n, edgeCount, dist, 100);
            
            double dijkstraTime = getAverageTime([&]() { dijkstra(testGraph, 0); });
            double deltaTime = getAverageTime([&]() { deltaStepping(testGraph, 0, DEFAULT_DELTA); });
            double parallelV1Time = getAverageTime([&]() { parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount); });
            double parallelV2Time = getAverageTime([&]() { parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount); });
            
            double sequentialSpeedup = dijkstraTime / deltaTime;
            double parallelV1Speedup = deltaTime / parallelV1Time;
            double parallelV2Speedup = deltaTime / parallelV2Time;
            
            results.push_back({
                dist,
                std::to_string(n),
                std::to_string(edgeCount),
                std::to_string(dijkstraTime),
                std::to_string(deltaTime),
                std::to_string(parallelV1Time),
                std::to_string(parallelV2Time),
                std::to_string(sequentialSpeedup),
                std::to_string(parallelV1Speedup),
                std::to_string(parallelV2Speedup)
            });
            
            std::cout << "  Size " << n << ": v1=" << std::fixed << std::setprecision(2) << parallelV1Speedup 
                      << "x, v2=" << parallelV2Speedup << "x speedup\n";
        }
    }
    
    writeResultsToFile("weight_distribution_analysis.csv", results);
    std::cout << "Weight distribution results saved\n\n";
}

void analyzeDeltaParameter() {
    std::cout << "Analyzing delta parameter sensitivity...\n";
    std::vector<std::vector<std::string>> results;
    results.push_back({"graph_type", "size", "delta", "optimal_delta", "delta_time", "parallel_v1_time", "parallel_v2_time", "v1_speedup", "v2_speedup"});
    
    std::vector<int> testSizes = {1000, 5000, 10000}; 
    std::vector<int> deltaValues = {1, 4, 6, 10, 15, 20, 30, 50};  
    std::vector<std::string> graphTypes = {"random", "grid", "scale_free"};
    
    int threadCount = std::thread::hardware_concurrency();
    
    for (const std::string& type : graphTypes) {
        std::cout << "Testing " << type << " graphs...\n";
        for (int n : testSizes) {
            Graph testGraph(1);
            if (type == "random") {
                testGraph = Graph::randomGraph(n, n * 2, 100);  
            } else if (type == "grid") {
                int side = (int)sqrt(n);
                testGraph = Graph::gridGraph(side, side, 50);
            } else if (type == "scale_free") {
                testGraph = Graph::scaleFreeGraph(n, 3, 100);
            }
            
            int optimal_delta = findDelta(testGraph);
            std::cout << "  Graph size: " << testGraph.size() << ", Optimal delta: " << optimal_delta << "\n";
            
            for (int delta : deltaValues) {
                double deltaTime = getAverageTime([&]() { deltaStepping(testGraph, 0, delta); });
                double parallelV1Time = getAverageTime([&]() { parallelDeltaStepping(testGraph, 0, delta, threadCount); });
                double parallelV2Time = getAverageTime([&]() { parallelDeltaStepping_v2(testGraph, 0, delta, threadCount); });
                
                double v1Speedup = deltaTime / parallelV1Time;
                double v2Speedup = deltaTime / parallelV2Time;
                
                results.push_back({
                    type,
                    std::to_string(testGraph.size()),
                    std::to_string(delta),
                    std::to_string(optimal_delta),
                    std::to_string(deltaTime),
                    std::to_string(parallelV1Time),
                    std::to_string(parallelV2Time),
                    std::to_string(v1Speedup),
                    std::to_string(v2Speedup)
                });
            }
            std::cout << "  Completed size " << testGraph.size() << "\n";
        }
    }
    
    writeResultsToFile("delta_parameter_analysis.csv", results);
    std::cout << "Delta parameter results saved\n\n";
}

void analyzeRealWorldGraphs() {
    std::cout << "Analyzing real-world graphs...\n";
    
    std::vector<std::vector<std::string>> results;
    results.push_back({"dataset", "vertices", "edges", "dijkstra_time", "delta_time", "parallel_v1_time", "parallel_v2_time", "seq_speedup", "par_v1_speedup", "par_v2_speedup", "correctness"});
    
    std::vector<std::string> graphFiles = {
        "network_graph/ca-GrQc.txt",
        "network_graph/ca-HepTh.txt", 
        "network_graph/ca-AstroPh.txt",
        "network_graph/ca-CondMat.txt",
        "network_graph/Slashdot0811.txt"
    };
    
    int threadCount = std::thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 4;
    
    for (const std::string& filename : graphFiles) {
        try {
            //first check correctness
            std::cout << "Processing " << filename << "...\n";
            Graph testGraph = Graph::loadSNAPGraph(filename);
            
            int edgeCount = 0;
            for (int i = 0; i < testGraph.size(); i++) {
                edgeCount += testGraph.getEdges(i).size();
            }
            edgeCount /= 2;
            std::cout << "  Running Dijkstra...\n";
            auto dijkstraResult = dijkstra(testGraph, 0);
            std::cout << "  Running sequential delta-stepping...\n";
            auto deltaResult = deltaStepping(testGraph, 0, DEFAULT_DELTA);
            std::cout << "  Running parallel v1...\n";
            auto parallelV1Result = parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount);
            std::cout << "  Running parallel v2...\n";
            auto parallelV2Result = parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount);
            bool isCorrect = true;
            std::string correctnessStatus = "PASS";
    
            if (dijkstraResult != deltaResult) {
                isCorrect = false;
                correctnessStatus = "FAIL:delta";
            } else if (dijkstraResult != parallelV2Result) {
                isCorrect = false;
                correctnessStatus = "FAIL:v2";
                for (int i = 0; i < dijkstraResult.size(); i++) {
                    if (dijkstraResult[i] != parallelV2Result[i]) {
                        std::cout << "  First mismatch at node " << i 
                                  << ": dijkstra=" << dijkstraResult[i] 
                                  << ", parallel_v2=" << parallelV2Result[i] << "\n";
                        break;
                    }
                }
            }
            if (!isCorrect) {
                std::cout << "  Correctness check failed: " << correctnessStatus << "\n";
                continue;
            }
            
            std::cout << "  ✓ Correctness is ok\n";
            
            // Now measure performance
            double dijkstraTime = getAverageTime([&]() { dijkstra(testGraph, 0); }, 3);
            double deltaTime = getAverageTime([&]() { deltaStepping(testGraph, 0, DEFAULT_DELTA); }, 3);
            double parallelV1Time = getAverageTime([&]() { parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threadCount); }, 3);
            double parallelV2Time = getAverageTime([&]() { parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threadCount); }, 3);
            
            double sequentialSpeedup = dijkstraTime / deltaTime;
            double parallelV1Speedup = deltaTime / parallelV1Time;
            double parallelV2Speedup = deltaTime / parallelV2Time;
            
            std::string datasetName = filename.substr(filename.find_last_of('/') + 1);
            
            results.push_back({
                datasetName,
                std::to_string(testGraph.size()),
                std::to_string(edgeCount),
                std::to_string(dijkstraTime),
                std::to_string(deltaTime),
                std::to_string(parallelV1Time),
                std::to_string(parallelV2Time),
                std::to_string(sequentialSpeedup),
                std::to_string(parallelV1Speedup),
                std::to_string(parallelV2Speedup),
                correctnessStatus
            });
            
            std::cout << "  " << testGraph.size() << " vertices: v1=" 
                      << std::fixed << std::setprecision(2) << parallelV1Speedup 
                      << "x, v2=" << parallelV2Speedup << "x speedup\n";
            
        } catch (const std::exception& e) {
            std::cout << "  Failed to load " << filename << ": " << e.what() << "\n";
        }
    }
    
    writeResultsToFile("real_world_analysis.csv", results);
    std::cout << "Real-world graph results saved\n\n";
}
void analyzeThreadScaling() {
    std::cout << "Analyzing thread scaling...\n";
    
    std::vector<std::vector<std::string>> results;
    results.push_back({"graph_type", "size", "threads", "delta_time", "parallel_v1_time", "parallel_v2_time", "v1_speedup", "v2_speedup"});
    
    int maxThreads = std::thread::hardware_concurrency();    
    std::vector<int> threadCounts = {1, 2, 4, maxThreads};
    
    std::vector<std::pair<std::string, int>> testConfigs = {
        {"scale_free_large", 25000},    
        {"dense_random", 15000},       
        {"grid_massive", 40000},        
        {"exponential_weights", 20000}  
    };
    
    for (const auto& [graphType, n] : testConfigs) {
        std::cout << "Testing " << graphType << " with " << n << " nodes...\n";
        
        Graph testGraph(1);
        
        if (graphType == "scale_free_large") {
            testGraph = Graph::scaleFreeGraph(n, 8, 100);
        }
        else if (graphType == "dense_random") {
            int edgeCount = n * 15;
            testGraph = Graph::randomGraph(n, edgeCount, 100);
        }
        else if (graphType == "grid_massive") {
            int side = (int)sqrt(n);
            testGraph = Graph::gridGraph(side, side, 50);
        }
        else if (graphType == "exponential_weights") {
            testGraph = Graph::randomGraphWithWeights(n, n * 8, "exponential", 100);
        }
        
        double baselineTime = getAverageTime([&]() { deltaStepping(testGraph, 0, DEFAULT_DELTA); });
        
        for (int threads : threadCounts) {
            double parallelV1Time = getAverageTime([&]() { 
                parallelDeltaStepping(testGraph, 0, DEFAULT_DELTA, threads); 
            });
            double parallelV2Time = getAverageTime([&]() { 
                parallelDeltaStepping_v2(testGraph, 0, DEFAULT_DELTA, threads); 
            });
            
            double v1Speedup = baselineTime / parallelV1Time;
            double v2Speedup = baselineTime / parallelV2Time;
            
            results.push_back({
                graphType,
                std::to_string(testGraph.size()),
                std::to_string(threads),
                std::to_string(baselineTime),
                std::to_string(parallelV1Time),
                std::to_string(parallelV2Time),
                std::to_string(v1Speedup),
                std::to_string(v2Speedup)
            });
            
            std::cout << "  " << threads << " threads: v1=" 
                      << std::fixed << std::setprecision(2) << v1Speedup 
                      << "x, v2=" << v2Speedup << "x speedup\n";
        }
        std::cout << "\n";
    }
    
    writeResultsToFile("thread_scaling_analysis.csv", results);
    std::cout << "Thread scaling results saved\n\n";
}

int main() {
    srand(static_cast<unsigned>(time(nullptr)));
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3);
    
    std::cout << "Delta-Stepping Algorithm Performance Analysis\n";
    std::cout << "============================================\n\n";
    
    if (!runCorrectnessTests()) {
        std::cerr << "Correctness tests failed - stopping analysis\n";
        return 1;
    }
    
    analyzeGraphTypes();
    // analyzeDensityEffects();
    // analyzeWeightDistributions();
    analyzeDeltaParameter();
    // analyzeRealWorldGraphs();
    // analyzeThreadScaling();
    
    std::cout << "Analysis complete! Generated files:\n";
    std::cout << "- graph_type_analysis.csv\n";
    std::cout << "- density_analysis.csv\n";
    std::cout << "- weight_distribution_analysis.csv\n";
    std::cout << "- delta_parameter_analysis.csv\n";
    std::cout << "- thread_scaling_analysis.csv\n";
    std::cout << "- real_world_analysis.csv\n";
    
    return 0;
}