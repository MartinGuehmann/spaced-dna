This is a fork of the *spaced words* approach to alignment-free dna comparision. Upstream is located at http://spaced.gobics.de/

# Installation

Regenerate the build system and compile:

	libtoolize --force
	aclocal
	autoheader
	automake --force-missing --add-missing
	autoconf
	./configure
	make

Execute with:

	./spaced [options] <sequence file>

The available options are:

	-h: print this help and exit
	-o <file>: output file name (default: DMat)
	-k <integer>: pattern weight (default 14)
	-l <integer>: pattern don't care positions (default 15)
	-n <integer>: number of patterns (default: 5)  
	-f <file>: use patterns in <file> instead of random patterns  
	-t <integer>: numer of threads (default: 25 threads)
	-r: don't consider the reverse complement
	-d EU | JS | EV: change distance type to Euclidean, Jensen-Shannon, evolutionary distance (default: EV) 

# References

*Spaced Words* is scientific software. If used in a scientific publication, please cite:

1. C.-A. Leimeister, M. Boden, S. Horwege, S. Lindner, B. Morgenstern (2014)
Fast alignment-free sequence comparison using spaced-word frequencies
Bioinformatics, DOI: 10.1093/bioinformatics/btu177 (http://bioinformatics.oxfordjournals.org/content/early/2014/04/03/bioinformatics.btu177)

1. S. Horwege, S. Linder, M. Boden, K. Hatje, M. Kollmar, C.-A. Leimeister, B. Morgenstern (2014)
Spaced words and kmacs: fast alignment-free sequence comparison based on inexact word matches
Nuc. Acids Research 42, W7-W11 (http://nar.oxfordjournals.org/content/42/W1/W7.abstract)

1. B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
Estimating evolutionary distances between genomic sequences from spaced-word matches
Algorithms for Molecular Biology 10, 5

# Contact

upstream contact:
chris.leimeister@stud.uni-goettingen.de
