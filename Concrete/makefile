# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20

# Targets and dependencies
main.out : main.o Matrix.o Matrix2.o MatrixN.o MatrixMxN.o Vector.o Vector2.o
	$(CXX) $(CXXFLAGS) -o main.out \
	Matrix.o Matrix2.o MatrixN.o MatrixMxN.o \
	Vector.o Vector2.o main.o

main.o : main.cpp Matrix2.o Vector2.o MatrixMxN.o
	$(CXX) $(CXXFLAGS) -c main.cpp

MatrixMxN.o : Matrix.o MatrixMxN.cpp
	$(CXX) $(CXXFLAGS) -c MatrixMxN.cpp

Matrix2.o : Matrix2.cpp MatrixN.o
	$(CXX) $(CXXFLAGS) -c Matrix2.cpp

MatrixN.o : MatrixN.cpp Matrix.o
	$(CXX) $(CXXFLAGS) -c MatrixN.cpp

Matrix.o : Matrix.hpp Matrix.cpp Vector.o
	$(CXX) $(CXXFLAGS) -c Matrix.cpp

Vector2.o : Vector.o Vector2.cpp
	$(CXX) $(CXXFLAGS) -c Vector2.cpp

Vector.o : Vector.hpp Vector.cpp
	$(CXX) $(CXXFLAGS) -c Vector.cpp

# Debug target
debug: CXXFLAGS += -g
debug: main.out

# Clean target
clear :
	rm -f *.o main.out

# Phony targets
.PHONY: clear debug
