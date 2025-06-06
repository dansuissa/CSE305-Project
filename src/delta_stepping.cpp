#include "delta_stepping.h"
#include <limits>
using std::vector;
static const int INF = std::numeric_limits<int>::max() / 2;

struct Request {
    int v;
    int x;
};

static void relax(int w, int x, vector<int>& dist, int delta, vector<vector<int>>& B) {
    if (x < dist[w]) {
        dist[w] = x;
        std::size_t j = static_cast<std::size_t>(x / delta);

        //difference from pseudocode: dynamic bucket array growth instead of cyclic reuse
        if (j >= B.size()) B.resize(j + 1);
        B[j].push_back(w);
    }
}

static vector<Request> findRequests(const Graph& g, const vector<int>& Vprime, int delta, const vector<int>& dist, bool lightPhase) {
    vector<Request> Rq;
    std::vector<Edge> light, heavy;
    for (int v : Vprime) {
        g.splitEdges(v, delta, light, heavy);
        const vector<Edge>& use = lightPhase ? light : heavy;
        for (const Edge& e : use) {
            int nd = dist[v] + e.weight;
            Rq.push_back({e.dest, nd});
        }
    }
    return Rq;
}

static void relaxRequests(const vector<Request>& Req, vector<int>& dist, int delta, vector<vector<int>>& B) {
    for (const Request& r : Req) relax(r.v, r.x, dist, delta, B);
}

vector<int> deltaStepping(const Graph& g, int source, int delta) {
    int n = g.size();
    int maxW = 0;
    for (int u = 0; u < n; ++u)
        for (const Edge& e : g.getEdges(u))
            if (e.weight > maxW) maxW = e.weight;
    std::size_t b = static_cast<std::size_t>(maxW / delta) + 1;
    vector<vector<int>> B(b);
    vector<int> dist(n, INF);
    relax(source, 0, dist, delta, B);
    std::size_t i = 0;
    while (true) {
        while (i < B.size() && B[i].empty()) ++i;
        if (i == B.size()) break;
        vector<int> R;
        while (!B[i].empty()) {
            vector<int> cur;
            cur.swap(B[i]);

            auto Req = findRequests(g, cur, delta, dist, true);
            R.insert(R.end(), cur.begin(), cur.end());
            relaxRequests(Req, dist, delta, B);
        }
        auto ReqH = findRequests(g, R, delta, dist, false);
        relaxRequests(ReqH, dist, delta, B);
    }
    return dist;
}
