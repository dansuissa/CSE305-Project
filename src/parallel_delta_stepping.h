#ifndef PARALLEL_DELTA_STEPPING_H
#define PARALLEL_DELTA_STEPPING_H
#include <vector>
#include "graph.h"
std::vector<int> parallelDeltaStepping(const Graph& g, int source, int delta = 0, int numThreads = 0);
#endif