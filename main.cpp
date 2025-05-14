#include "graph.h"
#include <iostream>

int main() {
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

    return 0;
}