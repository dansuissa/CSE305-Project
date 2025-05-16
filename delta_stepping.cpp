#include "delta_stepping.h"
#include <limits>
using std::vector;
static const int INF=std::numeric_limits<int>::max()/2;

vector<int> deltaStepping(const Graph& g,int src,int delta){
    int n=g.size();
    vector<int> dist(n,INF);
    vector<vector<int>> buckets(1);
    dist[src]=0;
    buckets[0].push_back(src);
    vector<Edge> light,heavy;
    size_t b=0;
    while(b<buckets.size()){
        if(buckets[b].empty()){++b;continue;}
        vector<int> old;
        while(!buckets[b].empty()){
            vector<int> cur;
            cur.swap(buckets[b]);
            for(int v:cur){
                old.push_back(v);
                g.splitEdges(v,delta,light,heavy);
                for(const Edge&e:light){
                    int w=e.dest;
                    int nd=dist[v]+e.weight;
                    if(nd<dist[w]){
                        dist[w]=nd;
                        size_t idx=nd/delta;
                        if(idx>=buckets.size())buckets.resize(idx+1);
                        buckets[idx].push_back(w);
                    }
                }
            }
        }
        for(int v:old){
            g.splitEdges(v,delta,light,heavy);
            for(const Edge&e:heavy){
                int w=e.dest;
                int nd=dist[v]+e.weight;
                if(nd<dist[w]){
                    dist[w]=nd;
                    size_t idx=nd/delta;
                    if(idx>=buckets.size())buckets.resize(idx+1);
                    buckets[idx].push_back(w);
                }
            }
        }
    }
    return dist;
}