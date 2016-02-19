/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * patternset object header
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
#ifndef PATTERNSET_H_
#define PATTERNSET_H_


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
#include "extkey.h"

class patternset {
public:
	patternset();
	patternset(char* pattern_file);
	patternset(int size, int *length, int weight);
	patternset(char* pattern_file, int size, int *length, int weight);
	~patternset();

	void ReInitPattern();

	void ChangeBits(int number);
	uint64_t ChangeBitPos(int pos, int pos_one, int pos_zero);
	bool UniqPattern(int number);

	std::vector<std::string> GetStringPattern();
	std::string GetString(int number);
	uint64_t GetPattern(int number);
	int GetWeight();
	int GetSize();
	int* GetLength();
	int GetLengthMean();
	bool GetImprove();
	bool GetUpdate();

	double GetValue(int number);
	void SetPattern(int number, uint64_t patt);
	void Print();
	void Silent();
	void RandPatLength();

protected:
	void Clear();
	void TestPattern();
	void VerifyConditions();

	std::vector<std::string> SplitString(std::string pattern, char* tokens);
	bool ValidatePatternsFormat(std::string pattern_form);
	bool ValidatePatternConditions();
	int PatternWeight(std::string pattern_wght);
	int* PatternLength(std::vector<std::string> pattern_length);

	void CreateLengths();
	void CreateRandomPattern();
	bool IsSetScore(int pattern);
	int SymbolRandPos(int number, char symb);
	int SymbolCalcPos(int number, char symb);
	std::vector<int> GetSymbol(int number, char symb);

	void ToString();
	uint64_t ToBit(std::string p);
	std::string BitString(uint64_t p);

	double MaxNumberPattern(int p_weight, int p_length);
	double Faculty(int value);
	void SecureMessage(std::string errmsg, int pos);

private:
	std::vector<uint64_t> pattern_set;
	std::vector<std::string> string_pat;
	std::vector<int> lengths;	
	int size;
	int *length;
	int length_mean;
	int weight;
	char* pattern_file = NULL;
	bool improve;
	bool update;
	bool silent;
	bool randpatleng;
};
#endif
