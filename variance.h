/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * variance object header
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
#ifndef VARIANCE_H_
#define VARIANCE_H_


#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "patternset.h"
#include "extkey.h"

class variance {
public:
	variance();
	variance(char* pattern_file);
	variance(int size, int *length, int weight);
	variance(char* pattern_file, int size, int *length, int weight, int l_hom, int l1, int l2, double p, double q, int H);
	~variance();

	void ReInitVariance();
	void RandPatLength();

	void Improve(int limit);
	void ImproveOC();

	std::vector<std::string> GetPattern();
	std::string GetPattern(int number);
	int GetWorstPat(int number);
	int GetWeight();
	int GetSize();
	int* GetLength();
	void PrintPattern();

	std::string GetFormat();
	double GetVariance();
	double GetNormVariance();
	double GetP();
	double GetQ();
	int GetLHom();
	int GetL1();
	int GetL2();
	int GetH();

	void LoopOpt(int n);
	void Quiet();
	void Silent();

protected:
	void Ctor(char* pattern_file, int size, int *length, int weight, int l_hom, int l1, int l2, double p, double q, int H);
	void InitVariance(char* pattern_file, int size, int *length, int weight);
	void InitVarMatrix();
	void InitVar();
	void SetPatOrder();
	void Clear();

	double Variance();
	void UpdateVariance(int pat);

	int ShiftPos(int number, int number2, int s);
	int WorstPattern(int number);
	void ResetPatternOrder();

	double Gauss();
	void SecureMessage(std::string errmsg);

private:
	std::vector<std::vector<double> > var_sum;
	patternset* pattern_set;
	std::vector<double> pattern_order;
	std::vector<extkey*> pattern_order_sort;
	std::string outvar;
	double variance_val;
	double p;
	double q;
	int size;
	int weight;
	int length_mean;
	int *length;
	int l_hom;
	int l1;
	int l2;
	int H;
	int loop_opt;
	bool oc;
	bool improve;
	bool quiet;
	bool silent;
	bool update; 
};
#endif
