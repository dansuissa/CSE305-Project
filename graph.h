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
    
    //for delta-stepping
    void splitEdges(int vertex, int delta, std::vector<Edge>& light, std::vector<Edge>& heavy) const;
    
    static Graph loadFile(const std::string& filename);
    
    void saveFile(const std::string& filename) const;
    
    static Graph randomGraph(int vertices, int edges, int maxWeight = 100);
    
    void print() const;
};

#endif


