/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * variance object file
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
#include "variance.h"


/*---Variables---------------------------------------------------------------*/



/*===Main-Part===============================================================*/
/*---Constructor-------------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues, pattern will be generated automatically.
 */
variance::variance() {
	int *length = new int[2];
	length[0] = 14;
	length[1] = 14;
	Ctor(NULL, 10, length, 8, 10000, 10000, 10000, 0.75, 0.25, 64);
}


/**
 * File constructor, sets values only from files; resets automatically if there are problems.
 */
variance::variance(char* pattern_file) {
	int *length = new int[2];
	length[0] = 14;
	length[1] = 14;
	Ctor(pattern_file, 10, length, 8, 10000, 10000, 10000, 0.75, 0.25, 64);
}


/**
 * Short constructor, sets some default vaulues, just pattern dimension is set.
 *
 * @param size		The amount of patterns; pattern number.
 *
 * @param length		The pattern length for each pattern of the pattern set.
 *
 * @param weigth		The weight (match positions; '1') for each pattern of the pattern set.
 */
variance::variance(int size, int *length, int weight) {
	Ctor(NULL, size, length, weight, 10000, 10000, 10000, 0.75, 0.25, 64);
}


/**
 * Long constructor, sets the values; resets automatically, if there are problems.
 *
 * @param pattern_file 	File, that may contains submitted pattern
 *
 * @param align_file 	File, that may contains an alignment file to estimate p, q, l_hom, li and lj
 *
 * @param size 		The amount of patterns; pattern number.
 *
 * @param length		The pattern length for each pattern of the pattern set.
 *
 * @param weigth		The weight (match positions; '1') for each pattern of the pattern set.
 *
 * @param l_hom		In theory, the amount of homologous positions of two sequences in an multiple alignment.
 *
 * @param l1		In theory, the first(represents each sequence i) sequence length of two observed sequences.
 *
 * @param l2		In theory, the second(represents each sequence j) sequence length of two observed sequences.
 *
 * @param p		The match probability ( = #matches / #l_hom)
 *
 * @param q		The background probability for each nucleotide A,C,G,T
 */
variance::variance(char* pattern_file, int size, int *length, int weight, int l_hom, int l1, int l2, double p, double q, int H) {
	Ctor(pattern_file, size, length, weight, l_hom, l1, l2, p, q, H);
}


/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
variance::~variance() {
	Clear();
}

void variance::Ctor(char* pattern_file, int size, int *length, int weight, int l_hom, int l1, int l2, double p, double q, int H){
	this->l_hom = l_hom;
	this->l1 = l1;
	this->l2 = l2;
	this->p = p;
	this->q = q;
	this->H = H;
	this->length = length;
	variance_val = 0;
	loop_opt = 1;
	quiet = false;
	silent = false;
	oc = false;
	update = false;
	outvar = "variance: ";
	InitVariance(pattern_file, size, length, weight);
}

/*---Init--------------------------------------------------------------------*/
/**
 * Complete Initializing of the Variance, which also creates the patternsets
 * 	calculates the first variance/oc and checks if each parameter is in its 
 * 	own correct size and dimension.
 */
void variance::InitVariance(char* pattern_file, int size, int *length, int weight){
	Clear();
	pattern_set = new patternset(pattern_file, size, length, weight);
	this->size = pattern_set->GetSize();
	this->weight = pattern_set->GetWeight();
	this->length = pattern_set->GetLength();
	length_mean = pattern_set->GetLengthMean();
	improve = pattern_set->GetImprove();
	update = pattern_set->GetUpdate();
	if (q > p || q >= 1 || p > 1 || q < 0 || p < 0) {						/*Incorrect values for variance have to be corrected*/
		SecureMessage("pq");
		p = 0.9;
		q = 0.25;
		update = true;
	}

	if (l_hom <= 0 || l1 <= 0 || l2 <= 0) {
		SecureMessage("length");
		l_hom = 10000;
		l1 = 10000;
		l2 = 10000;
		update = true;
	}
	InitVarMatrix();
	InitVar();
	if(update){
		Quiet();
	}
}

/**
 * Part Initialzinig of the Variance, it will not reset complete.
 * Only for ReInitializin a new pattternset with most values in
 * 	common.
 *  Public Access, so other programs can 'reload' it
 */
void variance::ReInitVariance(){
	for(int i = 0; i < size; i++){
		var_sum[i].clear();
	}
	var_sum.clear();
	pattern_order.clear();
	pattern_set->ReInitPattern();
	length_mean = pattern_set->GetLengthMean();
	InitVarMatrix();
	InitVar();
}


/**
 * Part of the destructor; could be used in other cases for 
 * 	reinitializing the variance/oc new.
 */
void  variance::Clear(){
	for(int i = 0; i < size; i++){
		var_sum[i].clear();
	}
	var_sum.clear();
	pattern_order.clear();
	delete pattern_set;
}

/**
 * Up to different patternsizes a specific matrix will be created.
 * Necessary for reinitliazing.
 */
void variance::InitVarMatrix() {
	std::vector<double> tmp;
	for (int i = 0; i < size; i++) {
		tmp.push_back(0.0);
	}
	for (int i = 0; i < size; i++) {
		var_sum.push_back(tmp);
	}
}

/**
 * First Initializing of the variance/oc where each value of the matrix has 
 * 	to be created, instead of updated variance, where it is not necessary to
 * 	recalculate each value/entry. Only of those, which alter due to pattern
 * 	modification.
 */
void variance::InitVar() {
	double var_hom, var_bac, hom, back, tmp;// ,hom_val, back_val;
	int shift, lower_length, upper_length;
	uint64_t oc_hom;

	oc_hom = 1;
	//hom_val = 0;
	//back_val = 0;

	for (int i = 0; i < size; i++) {								/*i and j represents Pi and Pj of the set of pattern*/
		lower_length = (int) log2(pattern_set->GetPattern(i));
		for (int j = i; j < size; j++) {
			upper_length = (int) log2(pattern_set->GetPattern(j)) + 1;

			if(i == j && !oc){
				upper_length = 1;
			}

			var_hom = 0.0;
			var_bac = 0.0;
			for (int s = -1 * lower_length; s < upper_length; s++) {			/*As in the formula, the shift goes from max shift left to max shift right*/
				shift = ShiftPos(i, j, s);						/*At least one position has to overlap*/
				if (!oc) {
					var_hom += (pow(p, shift) - pow(p, 2 * weight));		/*summation of the homologue first part*/
					var_bac += (pow(q, shift) - pow(q, 2 * weight));
				}
				else{	
					var_hom += (double) (oc_hom << (shift)); 				
				}
			}
			if(!oc){
				hom = (l_hom - length_mean + 1);
				hom *= var_hom;

				back = (l1 - length_mean + 1);
				back *= (l2 - length_mean);
				back *= var_bac;
			}
			else{
				hom = var_hom;
				back = 0.0; 
			}
		//	hom_val += hom;
		//	back_val += back;
			var_sum[i][j] = hom + back;							/*For each pair Pi and Pj this is the direct contribute to the complete variance*/
			var_sum[j][i] = var_sum[i][j];			
		}
	}
	variance_val = Variance();	

	pattern_order.clear();
	
	for(int i = 0; i < size; i++){
		tmp = 0.0;
		for(int j = 0; j < size; j++){
			tmp += var_sum[i][j];								/*Sets the contribute of each pattern in the beginning*/
		}
		pattern_order.push_back(tmp);
	}
	ResetPatternOrder();

	//std::cout << "homologue contribute: " << hom_val/variance_val << std::endl;
	//std::cout << "background contribute: " << back_val/variance_val << std::endl;
}


/*---Variance-----------------------------------------------------------------*/
/**
 * IF neccessary, it calculates the current variance/oc from the matrix.
 *
 * @return returns current variance/oc
 */
double variance::Variance() {
	double var;

	var = 0.0;

	for (int i = 0; i < size; i++) {
		for(int j = i; j < size; j++){
			var += var_sum[i][j];
		}
	}
	variance_val = var;
	return var;
}

/**
 * Instead of recalculating each value, only the matrix entries which
 * 	belong to the modified pattern have to be recalculated. 
 * 
 * @param pat	Patternindex of the modified pattern 
 */
void variance::UpdateVariance(int pat) {
	double var_hom, var_bac, hom, back;
	int shift, lower_length, upper_length;
	uint64_t oc_hom;	

	oc_hom = 1;

	for (int j = 0; j < size; j++) {
		var_hom = 0.0;
		var_bac = 0.0;

		lower_length = (int) log2(pattern_set->GetPattern(pat));
		upper_length = (int) log2(pattern_set->GetPattern(j)) + 1;		
		
		if(pat == j && !oc){
			upper_length = 1;
		}

		for (int s = -1 * lower_length; s < upper_length; s++) {				/*As in the formula, the shift goes from max shift left to max shift right*/
			shift = ShiftPos(pat, j, s);							/*At least one position has to overlap*/
			if (!oc) {
				var_hom += (pow(p, shift) - pow(p, 2 * weight));			/*summation of the homologue first part*/
				var_bac += (pow(q, shift) - pow(q, 2 * weight));
			}
			else{
				var_hom += (double) (oc_hom << (shift));			
			}
		}

		variance_val -= var_sum[pat][j];
		pattern_order[pat] -= var_sum[pat][j];
		pattern_order[j] -= var_sum[pat][j];
		
		if(!oc){
			hom = (l_hom - length_mean + 1);
			hom *= var_hom;
			back = (l1 - length_mean + 1);
			back *= (l2 - length_mean);
			back *= var_bac;
		}
		else{
			hom = var_hom;
			back = 0.0; 
		}

		var_sum[pat][j] = hom + back;								/*For each pair Pi and Pj this is the direct share of the complete variance...*/
		var_sum[j][pat] = var_sum[pat][j];

		variance_val += var_sum[pat][j];
		pattern_order[pat] += var_sum[pat][j];
		pattern_order[j] += var_sum[pat][j];
	}
}

/**
 * Shifts pattern and counts the number of common match positions
 *
 * @param p1 	Position of the first used pattern of the pattern set
 *
 * @param p2	Position of the second used pattern of the pattern set
 *		NOTE: possible is p1 = p2
 *
 * @param s	The shift of the second pattern, s < 0 := shift left pattern 2, s > := shift right pattern 2
 * 
 * @return 	Calculates and returns current variance
 */
int variance::ShiftPos(int number, int number2, int s) {
	int counter, maxi, maxa, maxb;
	uint64_t pata, patb, c;

	pata = pattern_set->GetPattern(number);
	patb = pattern_set->GetPattern(number2);

	counter = 0;
	c = 1;												/*LookUp Counter*/

	if (s < 0) {											/*Negative shift means, the upper pattern is shifted right*/
		s = 0 - s;

		pata = pata >> s;
	}
	else{												/*otherweise the lower pattern is shifted right*/
		patb = patb >> s;
	}

	maxa = (int) log2(pata)+1;
	maxb = (int) log2(patb)+1;
	maxi = std::min(maxa,maxb);
	for (int i = 0; i < maxi; i++) {
		if(((pata & patb) & c << i) != 0){
			counter++;
		}
	}

	if(!oc){
		counter = 2*weight - counter;								/*Allows us to calculate th n(p,p',s) value in a way...*/
	}												/*...with not so much space. We do not need to look for...*/
	return counter;											/*...each position containing at least one match*/
}

/**
 * Returns the position of the worst pattern, estimated by the maximum variancepart
 * 	for each pattern pair
 *
 * @param n		The pattern index of the chosen pattern.
 *
 * @return 		returns position worst matrix by max_value
 */
int variance::WorstPattern(int n){
	return pattern_order_sort[n]->GetPos();
}

void variance::ResetPatternOrder(){
	extkey* tmp;
	int end;

	if(pattern_order_sort.size() != 0){
		pattern_order_sort.clear();
	}

	if(length[0]==weight){
		end = size-1;
	}
	else{
		end = size;
	}

	for(int i = 0; i < end; i++){
		tmp = new extkey(i,pattern_order[i]);							/*Using an extended key to store each contribute and the ...*/
		pattern_order_sort.push_back(tmp);									/*pattern number.*/
	}

	std::sort(pattern_order_sort.rbegin(), pattern_order_sort.rend());	
}

/**
 * The improvement method
 *
 * @param limit		The number of patterns which have to be selected and mayby modified.
 */
void variance::Improve(int limit) {
	std::vector<uint64_t> vec;
	uint64_t pat_save;
	double var_save, tmp_var, best_by;
	int worst_pat, better_pattern, counter;
	//std::string str_out = "Output_variance";
	//std::ofstream fout;
	//fout.open(str_out);

	best_by = variance_val;

	if (!silent) {
		std::cout << "First " << outvar << "\t" << GetVariance() << std::endl;
		std::cout << "First norm_" << outvar << "\t" << GetNormVariance() << std::endl << std::endl;
	}
	//fout << 0 << "\t" << GetNormVariance() << std::endl;

	if(loop_opt > 1){
		pattern_set->Silent();
		for(int i = 0; i < size; i++){
			vec.push_back(pattern_set->GetPattern(i));
		}
	}


	for(int d = 0; d < loop_opt; d++){

		counter = 0;
		better_pattern = 0;
		worst_pat = 0;
		if(d > 0){
			if(!silent && !quiet){
				std::cout << "\n\n=================\nAnd do it again!" << std::endl;
			}
			ReInitVariance();
		}
		if (improve) {
			for (int i = 1; i <= limit; i++) {
				var_save = variance_val;
				worst_pat = GetWorstPat(counter%((int)pattern_order_sort.size()));
				//std::cout << counter << "\t" << worst_pat << std::endl;
				pat_save = pattern_set->GetPattern(worst_pat);
				pattern_set->ChangeBits(worst_pat);

				UpdateVariance(worst_pat);	
				tmp_var = variance_val;	

				if (tmp_var < var_save && pattern_set->UniqPattern(worst_pat)) {
					better_pattern++;
					var_save = variance_val;
					if (quiet && !silent) {
						std::cout << "\r*** BETTER PATTERN " << better_pattern << " ***";
						std::cout.flush();
					}
					else if (!quiet && !silent) {
						std::cout << "*** BETTER PATTERN " << better_pattern << " ***" << std::endl;
						std::cout << "Step " << i << " / " << limit << std::endl << "Patternset: \n";
						pattern_set->Print();
						std::cout << outvar << GetVariance() << std::endl;
						std::cout << "norm_" << outvar << GetNormVariance() << std::endl << std::endl;
					}
					else{
						//Do-Nothing
					}
					counter = 0;
					ResetPatternOrder();
					//fout << i << "\t" << GetNormVariance() << std::endl;
				}
				else{	
					pattern_set->SetPattern(worst_pat, pat_save);
					UpdateVariance(worst_pat);
					counter++;
				}
			}
			if(best_by > variance_val && loop_opt > 2){
				best_by = variance_val;
				for(int k = 0; k < size; k++){
					vec[k] = pattern_set->GetPattern(k);
				}
			}
			if (!silent && loop_opt < 2) {
				std::cout << "\n\nBest pattern:\n\n";
				pattern_set->Print();
				std::cout << "\nBest " << outvar << GetVariance() << std::endl;
				std::cout << "Best norm_"<< outvar << GetNormVariance() << std::endl;
			}
			//fout << limit << "\t" << GetNormVariance() << std::endl;
		}
	}
	//if(update){
		//SecureMessage("update");
	//}
	if(loop_opt > 1){
		for(int i = 0; i < size; i++){
			pattern_set->SetPattern(i, vec[i]);
		}
		variance_val = best_by;
		if(!silent){
			std::cout << "\n\nBest pattern:\n\n";
			pattern_set->Print();
			std::cout << "\nBest " << outvar << GetVariance() << std::endl;
			std::cout << "Best norm_"<< outvar << GetNormVariance() << std::endl;
		}
	}
	vec.clear();
	//fout.close();
}


/**
 * Turns the current calculation from variance to overlap complexity.
 * The matrix has to be recalculated complete.
 */
void variance::ImproveOC() {
	this->oc = true;
	outvar = "oc: ";
	InitVar();
}

/*---stuff-------------------------------------------------------------------*/
/**
 * Method calculate the number of all pattern combinations
 * Actually the gauss summation (n*(n+1)/2)
 *
 * @return the number of maximum pattern combinations
 */
double variance::Gauss() {
	return 0.5*size*size + 0.5*size;
}

/**
 * Changes to quiet output, only the number of better patterns will be printed 
 * 	automatically, and the best patternset by the variance/oc minimization.
 * Errors arr going to be printed!
 */
void variance::Quiet() {
	this->quiet = true;
}

/**
 * Changes to silent output, nothing will be printed automatically!
 * Errors are going to be printed!
 */
void variance::Silent() {
	this->quiet = true;
	this->silent = true;
	pattern_set->Silent();
}

/**
 * Resets the pattern lengths from evenly distributed lengths between maximum
 *	and minimum to random created pattern lengths.
 * One quarter of the set will have maximum length, the rest will have the
 *	length between maximum and minimum randomly distributed.
 */
void variance::RandPatLength(){
	pattern_set->RandPatLength();
	InitVar();
}

/**
 * Sets the number of times the improvement method should be executed again.
 * It will save the best variance/oc and its patternset.
 *
 * @param n		The number of iterations for the improvement method
 */
void variance::LoopOpt(int n){
	if(n == 0){
		n = 1;
	}
	loop_opt = n;
	
}

/**
 * Prints complete the current patternset
 */
void variance::PrintPattern(){
	pattern_set->Print();
}

/**
 * A Method to collect all errormessages. Just easier for programmer to change
 *  	the text or extend.
 *
 * @param errmsg
 *	Due to a few possible errormessages, this is the option, which has to be printed.
 *
 * @param pos
 *	The position of the incorrect patterns.
 *
 */
void variance::SecureMessage(std::string errmsg) {
	/*if (errmsg == "noimprove") {
		std::cerr << "Using your pattern conditions it is not sensible to improve your pattern, sorry!" << std::endl;
		std::cerr << "Deactivating improve mode\n" << std::endl;
		return;
	}
	if (errmsg == "update"){
		std::cerr << "\n\n--> IMPORTANT <--\nDue to some configuration errors, your submitted parameters have been updated!\n" << std::endl;
		return;
	}
	if (errmsg == "pq") {
		std::cerr << "Error while parsing your p and/or q value: \t0 < q <= p <= 1!" << std::endl;
		std::cerr << "Return to default values:\tp = 0.9 \tq=0.25\n" << std::endl;
		return;
	}
	if (errmsg == "length") {
		std::cerr << "Error while parsing your sequence length S: \t 0 < S!" << std::endl;
		std::cerr << "Return to default value: \tS = 10000\n" << std::endl;
		return;
	}
	if (errmsg == "wrongindex") {
		std::cerr << "ERROR! Pattern does not exist... do nothing\n" << std::endl;
		return;
	}
	if (errmsg == "speed"){
		std::cerr << "\n\n\nUPDATE! SpEED lenghts optimization is used, creating a new patternset!!!\n\n\n\n" << std::endl;
		return;
	}*/
}

/*---Set-&-GetFunc-----------------------------------------------------------*/

std::vector<std::string> variance::GetPattern(){
	return pattern_set->GetStringPattern();
}

std::string variance::GetPattern(int number){
	return pattern_set->GetString(number);
}

/**
 * Returns the current variance
 *
 * @return returns variance
 */
double variance::GetVariance() {
	return variance_val;
}

/**
 * Returns the current variance, normalized
 *
 * @return returns norm_variance
 */
double variance::GetNormVariance() {
	return variance_val / Gauss();
}


/**
 * Returns the position of the worst pattern, estimated by the maximum variancepart
 * 	for each pattern pair
 *
 * @return returns position worst matrix by max_value
 */
int variance::GetWorstPat(int number) {
	if(number >= size){
		return 0;
	}
	return WorstPattern(number);
}

/**
 * Returns the weight of each pattern, the match positions
 *
 * @return returns weight
 */
int variance::GetWeight() {
	return pattern_set->GetWeight();
}

/**
 * Returns the amount of patters; number of patterns
 *
 * @return returns size
 */
int variance::GetSize() {
	return size;
}

/**
 * Returns the length of each pattern
 *
 * @return returns length
 */
int* variance::GetLength() {
	length = pattern_set->GetLength();
	return length;
}

/**
 * Returns the match probability
 *
 * @return returns p value
 */
double variance::GetP() {
	return p;
}

/**
 * Returns the background probabillity, summation over all nucleotids
 *
 * @return returns q value
 */
double variance::GetQ() {
	return q;
}

/**
 * Returns the length of the homologous sequence pair
 *
 * @return returns homologous positions
 */
int variance::GetLHom() {
	return l_hom;
}

/**
 * Returns the length of the first observed sequence
 *
 * @return returns length sequence 1
 */
int variance::GetL1() {
	return l1;
}

/**
 * Returns the length of the second observed sequence
 *
 * @return 		returns length sequence 2
 */
int variance::GetL2() {
	return l2;
}

/**
 * Returns the current format, variance or overlap complexity
 *
 * @return 		the format string, either 'variance' or 'oc' 
 */
std::string variance::GetFormat() {
	return outvar;
}

/**
 * Returns the length of a random homolgue region on a dataset
 *
 * @return		The length of random homolgue region
 */
int variance::GetH(){
	return H;
}
