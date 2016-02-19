CPPFLAGS+=-Wall -Wextra -std=c++11
CXXFLAGS+=-O3 -fopenmp

%.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $^

spaced: spaced.o variance.o patternset.o extkey.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f *.o spaced
