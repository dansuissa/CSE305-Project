CXX = g++
CXXFLAGS = -std=c++17 -Wall -pthread

SRCS = main.cpp graph.cpp delta_stepping.cpp dijkstra.cpp \
       parallel_delta_stepping.cpp parallel_delta_stepping_v2.cpp

OBJS = $(SRCS:.cpp=.o)

main: $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)
clean:
	rm -f *.o main

all: main