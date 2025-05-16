CXX = g++
CXXFLAGS = -std=c++11 -Wall -pthread

all: main

main: main.o graph.o delta_stepping.o dijkstra.o
	$(CXX) -o main main.o graph.o delta_stepping.o dijkstra.o $(CXXFLAGS)

main.o: main.cpp graph.h delta_stepping.h dijkstra.h
	$(CXX) -c main.cpp $(CXXFLAGS)

graph.o: graph.cpp graph.h
	$(CXX) -c graph.cpp $(CXXFLAGS)

delta_stepping.o: delta_stepping.cpp delta_stepping.h graph.h
	$(CXX) -c delta_stepping.cpp $(CXXFLAGS)

dijkstra.o: dijkstra.cpp dijkstra.h graph.h
	$(CXX) -c dijkstra.cpp $(CXXFLAGS)

clean:
	rm -f *.o main
