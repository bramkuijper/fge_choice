CXX=g++
CXXFLAGS=-Wall -O3 -ggdb -std=c++14
EXE=solver.exe
CPP=solver.cpp solver_main.cpp
HPP=parameters.hpp solver.hpp

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

