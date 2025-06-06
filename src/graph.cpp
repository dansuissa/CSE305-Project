#include "graph.h"
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>

Edge::Edge(int dest, int weight) : dest(dest), weight(weight) {}

Graph::Graph(int n) : n(n) {
    edges.resize(n);
}

void Graph::addEdge(int from, int dest, int weight) {
    if (from >= 0 && from < n && dest >= 0 && dest < n) {
        edges[from].push_back(Edge(dest, weight));
    }
}

int Graph::size() const {
    return n;
}

const std::vector<Edge>& Graph::getEdges(int vertex) const {
    return edges[vertex];
}

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

Graph Graph::randomGraph(int vertices, int edges, int maxWeight) {
    Graph g(vertices);
    std::random_device rd;
    // std::mt19937 rng(rd());
    std::mt19937 rng(42); 
    std::uniform_int_distribution<> vertex_dist(0, vertices - 1);
    std::uniform_int_distribution<> weight_dist(1, maxWeight);
    
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

Graph Graph::gridGraph(int rows, int cols, int maxWeight) {
    Graph g(rows * cols);
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> weightDist(1, maxWeight);
    
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int node = r * cols + c;
            
            if (c < cols - 1) {
                int rightNode = r * cols + (c + 1);
                int weight = weightDist(gen);
                g.addEdge(node, rightNode, weight);
                g.addEdge(rightNode, node, weight);
            }
            
            if (r < rows - 1) {
                int bottomNode = (r + 1) * cols + c;
                int weight = weightDist(gen);
                g.addEdge(node, bottomNode, weight);
                g.addEdge(bottomNode, node, weight);
            }
        }
    }
    return g;
}

Graph Graph::scaleFreeGraph(int numNodes, int connectionsPerNode, int maxWeight) {
    Graph g(numNodes);
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> weightDist(1, maxWeight);
    
    int initialNodes = std::min(connectionsPerNode + 1, numNodes);
    for (int i = 1; i < initialNodes; i++) {
        int weight = weightDist(gen);
        g.addEdge(0, i, weight);
        g.addEdge(i, 0, weight);
    }
    
    std::vector<int> nodeDegrees(numNodes, 0);
    for (int i = 0; i < initialNodes; i++) {
        nodeDegrees[i] = 1;
    }
    
    for (int newNode = initialNodes; newNode < numNodes; newNode++) {
        std::vector<double> connectionProbs(newNode);
        double totalProb = 0;
        
        for (int existing = 0; existing < newNode; existing++) {
            connectionProbs[existing] = nodeDegrees[existing] + 1.0;
            totalProb += connectionProbs[existing];
        }
        
        std::uniform_real_distribution<> probDist(0.0, totalProb);
        std::set<int> chosenNodes;
        
        while (chosenNodes.size() < std::min(connectionsPerNode, newNode)) {
            double randomVal = probDist(gen);
            double cumulative = 0;
            
            for (int existing = 0; existing < newNode; existing++) {
                cumulative += connectionProbs[existing];
                if (randomVal <= cumulative && chosenNodes.find(existing) == chosenNodes.end()) {
                    chosenNodes.insert(existing);
                    break;
                }
            }
        }
        
        for (int target : chosenNodes) {
            int weight = weightDist(gen);
            g.addEdge(newNode, target, weight);
            g.addEdge(target, newNode, weight);
            nodeDegrees[newNode]++;
            nodeDegrees[target]++;
        }
    }
    return g;
}

Graph Graph::smallWorldGraph(int numNodes, int nearestNeighbors, double rewireProb, int maxWeight) {
    Graph g(numNodes);
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> weightDist(1, maxWeight);
    std::uniform_real_distribution<> probDist(0.0, 1.0);
    
    for (int i = 0; i < numNodes; i++) {
        for (int k = 1; k <= nearestNeighbors / 2; k++) {
            int neighbor = (i + k) % numNodes;
            int weight = weightDist(gen);
            g.addEdge(i, neighbor, weight);
            g.addEdge(neighbor, i, weight);
        }
    }
    
    for (int i = 0; i < numNodes; i++) {
        for (int k = 1; k <= nearestNeighbors / 2; k++) {
            if (probDist(gen) < rewireProb) {
                int newTarget;
                do {
                    newTarget = std::uniform_int_distribution<>(0, numNodes - 1)(gen);
                } while (newTarget == i);
                
                int weight = weightDist(gen);
                g.addEdge(i, newTarget, weight);
                g.addEdge(newTarget, i, weight);
            }
        }
    }
    return g;
}

Graph Graph::randomGraphWithWeights(int numNodes, int numEdges, const std::string& weightType, int maxWeight) {
    Graph g(numNodes);
    std::mt19937 gen(std::random_device{}());
    
    auto generateWeight = [&]() -> int {
        if (weightType == "uniform") {
            return std::uniform_int_distribution<>(1, maxWeight)(gen);
        } else if (weightType == "exponential") {
            std::exponential_distribution<> dist(2.0 / maxWeight);
            return std::max(1, std::min(maxWeight, (int)(dist(gen) * maxWeight)));
        } else if (weightType == "normal") {
            std::normal_distribution<> dist(maxWeight / 2.0, maxWeight / 6.0);
            return std::max(1, std::min(maxWeight, (int)dist(gen)));
        } else if (weightType == "bimodal") {
            if (std::uniform_real_distribution<>(0, 1)(gen) < 0.5) {
                return std::uniform_int_distribution<>(1, maxWeight / 4)(gen);
            } else {
                return std::uniform_int_distribution<>(3 * maxWeight / 4, maxWeight)(gen);
            }
        }
        return std::uniform_int_distribution<>(1, maxWeight)(gen);
    };
    
    for (int i = 1; i < numNodes; i++) {
        int parent = std::uniform_int_distribution<>(0, i - 1)(gen);
        int weight = generateWeight();
        g.addEdge(parent, i, weight);
        g.addEdge(i, parent, weight);
    }
    
    int edgesAdded = numNodes - 1;
    while (edgesAdded < numEdges) {
        int from = std::uniform_int_distribution<>(0, numNodes - 1)(gen);
        int to = std::uniform_int_distribution<>(0, numNodes - 1)(gen);
        
        if (from != to) {
            bool exists = false;
            for (const Edge& e : g.getEdges(from)) {
                if (e.dest == to) {
                    exists = true;
                    break;
                }
            }
            
            if (!exists) {
                int weight = generateWeight();
                g.addEdge(from, to, weight);
                edgesAdded++;
            }
        }
    }
    return g;
}

Graph Graph::loadSNAPGraph(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<int, int>> edgeList;
    std::set<int> nodeSet;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream stream(line);
        int from, to;
        if (stream >> from >> to) {
            edgeList.push_back({from, to});
            nodeSet.insert(from);
            nodeSet.insert(to);
        }
    }
    file.close();
    
    std::map<int, int> nodeMapping;
    int newId = 0;
    for (int originalId : nodeSet) {
        nodeMapping[originalId] = newId++;
    }
    
    Graph g(nodeSet.size());
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> weightDist(1, 100);
    
    for (const auto& edge : edgeList) {
        int mappedFrom = nodeMapping[edge.first];
        int mappedTo = nodeMapping[edge.second];
        int weight = weightDist(gen);
        g.addEdge(mappedFrom, mappedTo, weight);
        g.addEdge(mappedTo, mappedFrom, weight);
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