#include "dijkstra.h"
#include <queue>
#include <limits>
using std::vector;
static const int INF=std::numeric_limits<int>::max()/2;

struct Node{
    int v,d;
    bool operator>(const Node& o)const{return d>o.d;}
};

vector<int> dijkstra(const Graph& g,int src){
    int n=g.size();
    vector<int> dist(n,INF);
    std::priority_queue<Node,std::vector<Node>,std::greater<Node>> pq;
    dist[src]=0;
    pq.push({src,0});
    while(!pq.empty()){
        Node cur=pq.top();pq.pop();
        if(cur.d!=dist[cur.v])continue;
        for(const Edge&e:g.getEdges(cur.v)){
            int w=e.dest;
            int nd=cur.d+e.weight;
            if(nd<dist[w]){
                dist[w]=nd;
                pq.push({w,nd});
            }
        }
    }
    return dist;
}