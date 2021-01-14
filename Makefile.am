bin_PROGRAMS = spaced
spaced_SOURCES = src/patternset.cpp src/patternset.h src/sort.h src/spaced.cc src/variance.cpp src/variance.h
spaced_CPPFLAGS = $(OPENMP_CPPFLAGS) -Wall -Wextra -std=c++11 -Wno-char-subscripts
spaced_CXXFLAGS = $(OPENMP_CXXFLAGS)

man_MANS=spaced.1

dist_noinst_DATA=README.md
