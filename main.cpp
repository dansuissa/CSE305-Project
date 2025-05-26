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

using Clock = std::chrono::steady_clock;
constexpr int DELTA_DEFAULT = 3;

template <typename F>
double measureMs(F &&f)
{
    auto t0 = Clock::now();
    f();
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

void showDistances(const std::vector<int> &d)
{
    for (std::size_t i = 0; i < d.size(); ++i)
        std::cout << i << ": " << d[i] << '\n';
}

// sanity check 
void testSmallGraph()
{
    std::cout << "sanity test:\n";
    Graph g(5);
    g.addEdge(0, 1, 10);
    g.addEdge(0, 2, 5);
    g.addEdge(1, 2, 2);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 9);
    g.addEdge(2, 4, 2);
    g.addEdge(3, 4, 4);
    std::cout << "graph:\n";
    g.print();
    double tDijkstra{0.0}, tDelta{0.0}, tPar{0.0};
    std::vector<int> dj, ds, par;
    int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    tDijkstra = measureMs([&]() { dj = dijkstra(g, 0); });
    tDelta    = measureMs([&]() { ds = deltaStepping(g, 0, DELTA_DEFAULT); });
    // tPar      = measureMs([&]() { par = parallelDeltaStepping(g, 0, DELTA_DEFAULT, threads); });
    tPar      = measureMs([&]() { par = parallelDeltaStepping_v2(g, 0, DELTA_DEFAULT, threads); });

    std::cout << "\nDijkstra (" << tDijkstra << " ms)\n";
    showDistances(dj);

    std::cout << "\n delta‑stepping (delta=" << DELTA_DEFAULT << ", " << tDelta << " ms)\n";
    showDistances(ds);

    std::cout << "\nparallel delta‑stepping (" << tPar << " ms)\n";
    showDistances(par);

    if (dj != ds || dj != par)
        std::cerr << "algorithm error\n";
    std::cout << '\n';
}


//benchmark

void performanceBenchmarks()
{
    std::cout << "scaling experiment:\n";
    std::cout << "  n    m      Dijkstra    delta‑step    parallel\n"
                 "---------------------------------------------\n";
    const int sizes[] = {50, 100, 200, 500};
    int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    for (int n : sizes)
    {
        int m = n * 3;
        Graph g = Graph::randomGraph(n, m);

        double tDijkstra = measureMs([&]() { dijkstra(g, 0); });
        double tDelta    = measureMs([&]() { deltaStepping(g, 0, DELTA_DEFAULT); });
        // double tPar      = measureMs([&]() { parallelDeltaStepping(g, 0, DELTA_DEFAULT, threads); });
        double tPar      = measureMs([&]() { parallelDeltaStepping_v2(g, 0, DELTA_DEFAULT, threads); });

        std::cout << std::setw(4) << n << ' '
                  << std::setw(5) << m << "   "
                  << std::setw(10) << std::fixed << std::setprecision(3) << tDijkstra << "   "
                  << std::setw(10) << tDelta << "   "
                  << std::setw(10) << tPar << '\n';
    }
    std::cout << '\n';
}


// cross validation
void CrossValidation()
{
    std::cout << "random validation:\n";
    bool ok = true;
    int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    for (int t = 0; t < 30 && ok; ++t)
    {
        int v = 20 + std::rand() % 50;
        int e = v * (2 + std::rand() % 3);
        Graph g = Graph::randomGraph(v, e);
        auto a = deltaStepping(g, 0, DELTA_DEFAULT);
        auto b = dijkstra(g, 0);
        // auto c = parallelDeltaStepping(g, 0, DELTA_DEFAULT, threads);
        auto c = parallelDeltaStepping_v2(g, 0, DELTA_DEFAULT, threads);
        if (a != b || a != c)
        {
            ok = false;
            std::cerr << "mismatch on graph " << t
                      << " (" << v << " v, " << e << " e)\n";
        }
    }
    if (ok)
        std::cout << "the parallel version is valid";
    std::cout << '\n';
}

void largeGraphBenchmark()
{
    std::cout << "Large graph performance test:\n";
    const int V = 50000;
    const int E = 5000000;
    std::cout << "random graph with " << V << " nodes and " << E << " edges...\n";

    Graph g = Graph::randomGraph(V, E);

    int hw = std::thread::hardware_concurrency();
    int threads = hw ? hw : 1;
    int optimal_delta = findDelta(g);

    double tSeq = measureMs([&] { deltaStepping(g, 0, DELTA_DEFAULT); });
    // double tPar = measureMs([&] { parallelDeltaStepping(g, 0, DELTA_DEFAULT, threads); });
    // double tAuto = measureMs([&] { parallelDeltaStepping(g, 0, 0, threads); }); // auto delta

    double tPar = measureMs([&] { parallelDeltaStepping_v2(g, 0, DELTA_DEFAULT, threads); });
    double tAuto = measureMs([&] { parallelDeltaStepping_v2(g, 0, optimal_delta, threads); }); // auto delta

    std::cout << "Sequential Delta-stepping: " << tSeq << " ms\n";
    std::cout << "Parallel Delta-stepping (" << threads << " threads, delta="
              << DELTA_DEFAULT << "): " << tPar << " ms\n";
    std::cout << "Auto-tuned Parallel Delta-stepping (delta= " << threads << ", threads="
              << threads << "): " << tAuto << " ms\n";

    std::cout << "Parallel Speedup: " << tSeq / tPar << "x\n";
    std::cout << "Auto-tuned Speedup: " << tSeq / tAuto << "x\n\n";
}

int main()
{
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3);
    testSmallGraph();
    performanceBenchmarks();
    CrossValidation();
    largeGraphBenchmark();
    return 0;
}
