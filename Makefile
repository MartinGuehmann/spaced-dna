CPPFLAGS=-Wall -Wextra -std=c++11
CXXFLAGS=-O3 -fopenmp


spaced: spaced.cc variance.cpp patternset.cpp extkey.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

clean:
	rm -f *.o spaced
