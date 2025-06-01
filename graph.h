#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <limits>

struct Edge {
    int dest;
    int weight;
    Edge(int dest, int w);
};

class Graph {
private:
    int n;
    std::vector<std::vector<Edge>> edges;

public:
    Graph(int n);
    
    void addEdge(int from, int to, int weight);
    int size() const;
    const std::vector<Edge>& getEdges(int vertex) const;
    void splitEdges(int vertex, int delta, std::vector<Edge>& light, std::vector<Edge>& heavy) const;
    static Graph loadFile(const std::string& filename);
    void saveFile(const std::string& filename) const;
    static Graph randomGraph(int vertices, int edges, int maxWeight = 100);
    static Graph gridGraph(int rows, int cols, int maxWeight = 10);
    static Graph scaleFreeGraph(int numNodes, int connectionsPerNode = 3, int maxWeight = 100);
    static Graph smallWorldGraph(int numNodes, int nearestNeighbors = 4, double rewireProb = 0.3, int maxWeight = 50);
    static Graph randomGraphWithWeights(int numNodes, int numEdges, const std::string& weightType, int maxWeight = 100);
    static Graph loadSNAPGraph(const std::string& filename);
    
    void print() const;
};

#endif