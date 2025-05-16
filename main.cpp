#include "graph.h"
#include "delta_stepping.h"
#include "dijkstra.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

int main() {
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    Graph g(5);
    g.addEdge(0, 1, 10);
    g.addEdge(0, 2, 5);
    g.addEdge(1, 2, 2);
    g.addEdge(1, 3, 1);
    g.addEdge(2, 3, 9);
    g.addEdge(2, 4, 2);
    g.addEdge(3, 4, 4);
    
    std::cout << "Test graph:\n";
    g.print();
    
    std::vector<Edge> light, heavy;
    int delta = 3;
    g.splitEdges(2, delta, light, heavy);
    std::cout << "\nTest splitting edges from vertex 2 with delta = " << delta << ":\n";
    std::cout << "Light edges: ";
    for (const Edge& e : light) {
        std::cout << e.dest << "(" << e.weight << ") ";
    }
    
    std::cout << "\nHeavy edges: ";
    for (const Edge& e : heavy) {
        std::cout << e.dest << "(" << e.weight << ") ";
    }
    std::cout << "\n";

    std::cout << "\nRandom graph test:\n";
    Graph random = Graph::randomGraph(8, 15);
    random.print();

    random.saveFile("test_graph.txt");
    Graph loaded = Graph::loadFile("test_graph.txt");
    loaded.print();

    int delta_ds = 3;
    auto dist_ds = deltaStepping(g, 0, delta_ds);
    std::cout << "\nDelta-stepping distances from 0 (delta = " << delta_ds << "):\n";
    for (size_t i = 0; i < dist_ds.size(); ++i) std::cout << i << ": " << dist_ds[i] << "\n";

    auto dist_dj = dijkstra(g, 0);
    std::cout << "\nDijkstra distances from 0:\n";
    for (size_t i = 0; i < dist_dj.size(); ++i) std::cout << i << ": " << dist_dj[i] << "\n";

    bool ok = true;
    for (int t = 0; t < 500 && ok; ++t) {
        int v = 20 + std::rand() % 81;
        int e = v * (2 + std::rand() % 4);      
        int max_w = 100;
        Graph rg = Graph::randomGraph(v, e, max_w);
        auto a = deltaStepping(rg, 0, 3);
        auto b = dijkstra(rg, 0);
        if (a != b) {
            std::cout << "\nMismatch on random test " << t << "\n";
            ok = false;
        }
    }
    if (ok) std::cout << "\nRandom-graph cross-check passed (500 cases)\n";

    return 0;
}
