CXX=g++
CXXFLAGS=-Wall -O3 -ggdb -std=c++17
IFLAGS=/usr/include/eigen3

solver.exe : solver_main.cpp solver.cpp parameters.hpp
	$(CXX) $(CXXFLAGS) -I $(IFLAGS) -o solver.exe solver_main.cpp solver.cpp
	
