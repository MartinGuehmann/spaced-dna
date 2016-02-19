/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * extended key object header
 *
 * For theory please have a look at:
 *
 * B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
 * Estimating evolutionary distances between genomic sequences from spaced-word matches
 * Algorithms for Molecular Biology 10, 5. (http://www.almob.org/content/10/1/5/abstract)
 *
 *
 * @author: Lars Hahn - 26.10.2015, Georg-August-Universitaet Goettingen
 * @version: 1.0.2 11/2015
 */
#ifndef EXTKEY_H_
#define EXTKEY_H_

class extkey{
	public:
		extkey();
		extkey(int pos, double value);

		int GetPos();
		double GetValue();
		void SetValue(double value);
	
		bool operator < (extkey const& p) const;
		bool operator <= (extkey const& p) const;
		bool operator > (extkey const& p) const;
		bool operator >= (extkey const& p) const;
		bool operator == (extkey const& p) const;
		bool operator != (extkey const& p) const;
		bool operator || (extkey const& p) const;

	private:
		int pos;
		double value;
};
#endif
