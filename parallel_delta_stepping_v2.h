#ifndef PARALLEL_DELTA_STEPPING_V2_H
#define PARALLEL_DELTA_STEPPING_V2_H

#include <vector>
#include "graph.h"


int findDelta(const Graph& g);
int optimalNumberOfThreads(const Graph& g, int maxThreads);
std::vector<int> parallelDeltaStepping_v2(const Graph& g, int source, int delta = -1, int numThreads = -1);

#endif // PARALLEL_DELTA_STEPPING_V2_H