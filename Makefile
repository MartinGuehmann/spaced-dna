
all: 
		g++ -O3 -fopenmp -std=c++11 src/spaced.cc  src/variance.cpp src/patternset.cpp -o spaced
 
