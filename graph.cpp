#include "graph.h"
#include <iostream>
#include <fstream>
#include <random>

Edge::Edge (int dest, int weight) : dest (dest), weight(weight) {}

Graph::Graph (int n) : n(n) {edges.resize(n);}

void Graph::addEdge(int from, int dest, int weight) {
    if (from >= 0 && from < n && dest >= 0 && dest < n) {
        edges[from].push_back(Edge(dest, weight));
    }
}

int Graph::size() const {return n;}

const std::vector<Edge>& Graph::getEdges(int vertex) const {return edges[vertex];}

void Graph::splitEdges(int vertex, int delta, std::vector<Edge>& light, std::vector<Edge>& heavy) const {
    light.clear();
    heavy.clear();
    for (const Edge& e : edges[vertex]) {
        if (e.weight <= delta) {
            light.push_back(e);
        } else {
            heavy.push_back(e);
        }
    }
}

Graph Graph::loadFile(const std::string& filename) {
    std::ifstream file(filename);
    int vcount, ecount;
    file >> vcount >> ecount;
    Graph g(vcount);
    for (int i = 0; i < ecount; i++) {
        int from, to, w;
        file >> from >> to >> w;
        g.addEdge(from, to, w);
    }
    return g;
}

void Graph::saveFile(const std::string& filename) const {
    std::ofstream file(filename);
    
    int total = 0;
    for (const auto& list : edges) {
        total += list.size();
    }
    
    file << n << " " << total << "\n";
    for (int i = 0; i < n; i++) {
        for (const Edge& e : edges[i]) {
            file << i << " " << e.dest << " " << e.weight << "\n";
        }
    }
}

Graph Graph::randomGraph(int vertices, int edges, int max_w) {
    Graph g(vertices);
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<> vertex_dist(0, vertices - 1);
    std::uniform_int_distribution<> weight_dist(1, max_w);
    for (int i = 1; i < vertices; i++) {
        int from = vertex_dist(rng) % i; 
        g.addEdge(from, i, weight_dist(rng));
    }
    int added = vertices - 1; 
    while (added < edges) {
        int from = vertex_dist(rng);
        int dest = vertex_dist(rng);
        if (from != dest) {
            bool exists = false;
            for (const Edge& e : g.edges[from]) {
                if (e.dest == dest) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                g.addEdge(from, dest, weight_dist(rng));
                added++;
            }
        }
    }
    return g;
}

void Graph::print() const {
    std::cout << "Graph with " << n << " vertices:\n";
    
    for (int i = 0; i < n; i++) {
        std::cout << i << " -> ";
        bool first = true;
        for (const Edge& e : edges[i]) {
            if (!first) std::cout << ", ";
            std::cout << e.dest << "(" << e.weight << ")";
            first = false;
        }
        
        std::cout << "\n";
    }
}