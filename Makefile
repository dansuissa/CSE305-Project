CXX = g++
CXXFLAGS = -std=c++11 -Wall -pthread

all: main

main: main.o graph.o
	$(CXX) -o main main.o graph.o $(CXXFLAGS)

main.o: main.cpp graph.h
	$(CXX) -c main.cpp $(CXXFLAGS)

graph.o: graph.cpp graph.h
	$(CXX) -c graph.cpp $(CXXFLAGS)

clean:
	rm -f *.o main