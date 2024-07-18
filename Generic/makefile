# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20

# Targets and dependencies
main.out : main.o
	$(CXX) $(CXXFLAGS) -o main.out main.o

main.o : main.cpp matrix.hpp vector.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

# Debug target
debug: CXXFLAGS += -g
debug: main.out

# Clean target
clear :
	rm -f *.o main.out

# Phony targets
.PHONY: clear debug
