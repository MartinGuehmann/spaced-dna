/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * patternset object file
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
#include "patternset.h"


/*---Variables---------------------------------------------------------------*/

std::default_random_engine generator(std::random_device{}());

/*===Main-Part===============================================================*/
/*---Constructor-------------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues, pattern will be generated automatically.
 */
patternset::patternset() {
	this->size = 10;
	this->length = new int[2];
	this->length[0] = 14;
	this->length[1] = 14;
	this->weight = 8;
	patternset(NULL, size, length, weight);
}


/**
 * File constructor, sets values only from files; resets automatically if there are problems.
 */
patternset::patternset(char* pattern_file) {
	this->pattern_file = pattern_file;
	this->length = new int[2];
	this->length[0] = 14;
	this->length[1] = 14;
	patternset(pattern_file, 10, length, 8);
}

/**
 * Short constructor, sets some default vaulues, just pattern dimension is set.
 *
 * @param size
 * 		The amount of patterns; pattern number.
 *
 * @param length
 *		The pattern length for each pattern of the pattern set.
 *
 * @param weigth
 * 		The weight (match positions; '1') for each pattern of the pattern set.
 */
patternset::patternset(int size, int *length, int weight) {
	patternset(NULL, size, length, weight);
}


/**
 * Long constructor, sets the values; resets automatically, if there are problems.
 *
 * @param pattern_file		File, that may contains submitted pattern
 *
 * @param align_file		File, that may contains an alignment file to estimate p, q, l_hom, li and lj
 *
 * @param size			The amount of patterns; pattern number.
 *
 * @param length			The pattern length for each pattern of the pattern set.
 *
 * @param weigth			The weight (match positions; '1') for each pattern of the pattern set.
 *
 * @param l_hom			In theory, the amount of homologous positions of two sequences in an multiple alignment.
 *
 * @param l1			In theory, the first(represents each sequence i) sequence length of two observed sequences.
 *
 * @param l2			In theory, the second(represents each sequence j) sequence length of two observed sequences.
 *
 * @param p			The match probability ( = #matches / #l_hom)
 *
 * @param q			The background probability for each nucleotide A,C,G,T
 */
patternset::patternset(char* pattern_file, int size, int *length, int weight) {
	this->pattern_file = pattern_file;
	this->size = size;
	this->length = length;
	this->weight = weight;
	update = false;
	improve = true;
	silent = false;
	randpatleng = false;
	ReInitPattern();
}

/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
patternset::~patternset() {
	delete[] length;
	Clear();
}

/*---Init--------------------------------------------------------------------*/

/**
 * Creates for the submitted or default values a set of pattern and calculates the first variance.
 * If possible estimates p, q, all sequences lengths and all combinations of homologous sequence positions
 *
 * Also reset Pattern if needed, not necessary to create new object
 */
void patternset::ReInitPattern() {			/*====Main-TODO=====*/
	uint64_t tmp;
	int lgth, lgt;
	Clear();
	
	lgth = 0;

	if (pattern_file != NULL) {
		TestPattern();
	}	
	if(string_pat.size()== 0){
		VerifyConditions();
		CreateRandomPattern();
	}
	else{
		lengths.clear();
		for(int i = 0; i < size; i++){
			tmp = ToBit(string_pat[i]);
			pattern_set.push_back(tmp);
			lgt = (int)(log2(tmp))+1;
			lengths.push_back(lgt);
		}
		length[0] = lengths[0];
		length[1] = length[0];
		for(int i = 1; i < size; i++){
			if(length[0] > lengths[i]){
				length[0] = lengths[i];
			}
			if(length[1] < lengths[i]){
				length[1] = lengths[i];
			}
		}
	}
	for(int i = 0; i < size; i++){
		lgth += lengths[i];
	}
	length_mean = lgth/size;	

	return;
}
/**
 * Due to some cases, it is possible necessary to
 * reset some pattern parameters.
 *
 * This is not a complete delete, so therefore it
 * it is only a part of the destructor.
 */
void patternset::Clear(){
	lengths.clear();
	pattern_set.clear();
	string_pat.clear();
}


/*---Functions---------------------------------------------------------------*/
/**
 * If there is an submitted pattern_file, it will check, if it is containing
 *	a possible pattern, and if this pattern is in the correct format.
 * Pattern can be parsed in a different formats
 */
void patternset::TestPattern() {
	std::ifstream patternfile;
	std::vector<std::string> pattern_tmp;
	std::string tmp;
	char tokens[4] = { '.',' ',',',';' };							/*These tokens are allowed to seperate patterns*/
	size_t f_size;
	bool start;

	start = false;

	patternfile.open(pattern_file);
	patternfile.seekg(0, std::ios::end);
	f_size = (size_t)patternfile.tellg();
	if (!patternfile) {
		SecureMessage("file", -1);
	}
	else if (f_size == 0) {									/*[FILE].eof() does not recognize empty files...*/
		SecureMessage("empty", -1);
	}
	else {
		patternfile.close();								/*..therefore it has also to be closed and opened -.-** */
		patternfile.open(pattern_file);
		if(!silent){
			std::cout << "Reading pattern from submitted patternfile ...\n" << std::endl;
		}
		patternfile >> tmp;
		while (!patternfile.eof()) {
			if (!ValidatePatternsFormat(tmp)) {
				patternfile >> tmp;						/*Ignoring each incorrect pattern. It is easier to calculate with all the rests*/
			}
			else {
				pattern_tmp = SplitString(tmp, tokens);
				for (unsigned int i = 0; i < pattern_tmp.size(); i++) {
					string_pat.push_back(pattern_tmp[i]);			/*Each pattern needs to be saved in its own std::string for comparison*/
				}
				patternfile >> tmp;
			}
		}
		if (size > 0) {									/*For an empty set we do not habe to validate*/
			start = ValidatePatternConditions();
		}
		if (!start) {									/*Some conditions, like changed weight is too much much amount of work*/
			SecureMessage("pattern", -1);						/*Therefore we just use our default values*/
			size = 10;
			weight = 8;
			length = new int[2];
			length[0] = 14;
			length[1] = 14;
			string_pat.clear();							/*Do not forget to reset, or the pattern set will not replaces, just increased*/
		}
		else {
			size = (int )string_pat.size();
			length = PatternLength(string_pat);
			weight = PatternWeight(string_pat[0]);
			if(!silent){
				for (int i = 0; i < size - 1; i++) {
					std::cout << string_pat[i] << "\n";
				}
				std::cout << string_pat[size - 1] << std::endl;
				std::cout << "\n... Done!\n" << std::endl;
			}
		}
	}
	patternfile.close();
}


/**
 * It will Verify the pattern conditions: is the weight correct,
 * 	not negative values, weight not above length and so on.. 
 */
void patternset::VerifyConditions(){
	std::vector<int> leng_new;
	int tmp, diff, leng_old;

	if(size <= 0){
		SecureMessage("size", -1);
		size = 10;
		update = true;
	}
	if(length[1] < length[0]){
		tmp = length[1];
		length[1] = length[0];
		length[0] = tmp;
	}	
	if(weight < 2){
		SecureMessage("weight_pat", -1);
		weight = 8;
		update = true;
	}

	if(length[0] < 3 || length[1] < 3){
		SecureMessage("length",-1);
		diff = length[1] - length[0];
		length[0] = weight*2;
		length[1] = length[0]+diff;
		update = true;
	}

	diff = length[1]-length[0];
	leng_old = length[0];
	if(length[1] > 63){								/*max length of 64-Bit Integer*/
		SecureMessage("length",-1);
		if(diff >= 63-weight){
			diff = 63-weight-1;
		}
		length[1] = 63;
		length[0] = length[1]-diff;
		update = true;
	}
	if(length[0] < weight){
		SecureMessage("length",-1);
		if(length[0] != leng_old){
			length[0] = leng_old;
		}
		else{
			length[0] = weight+2;
		}
		if(length[0] > length[1]){
			length[0] = length[1];
		}
		update = true;
	}
	
	if(length[1] < weight || length[0] < weight || length[0] > 63 || length[1] > 63){
		SecureMessage("nkl", -1);					/*...if everything is fucked up and rescue does not work...*/
		length[0] = 14;
		length[1] = 14;
		weight = 7;
		update = true;
	}

	if(length[0] == weight && length[1] == weight){
		if(size > 1){
			SecureMessage("max_number_pattern",-1);
			size = 1;
			update = true;
		}
		if(improve){
			improve = false;
			SecureMessage("noimprove",-1);
			update = true;
		}
	}

	CreateLengths();

	tmp = 1;
	diff = 0;
	if(lengths[0] == weight && size > 1){
		diff++;
		if(lengths[1] == lengths[0]){
			lengths[1]++;
		}
	}

	for(int i = 1+diff; i < size; i++){
		if(lengths[i] < lengths[i-1]){
			lengths[i] = lengths[i-1];
		}
		if(lengths[i] == lengths[i-1]){
			tmp++;
			if(tmp > MaxNumberPattern(weight-2, lengths[i]-2)){
				lengths[i]++;
				if(lengths[i] > length[1]){
					lengths[i] = length[1];
					size = i;
					if(size <= 0){
						size = 1;
					}
				}
				tmp = 1;
			}
		}
		else{
			tmp = 1;
		}
	}


	if(size != (int)lengths.size()){
		SecureMessage("max_number_pattern",-1);
		for(int i = 0; i < size; i++){
			leng_new.push_back(lengths[i]);
		}
		lengths.clear();
		lengths = leng_new;
		update=true;
	}
	return;
}



/*---Functions---------------------------------------------------------------*/
/*---Func-Create-------------------------------------------------------------*/
/**
 * Splits a string, read by a pattern file. Mayby each pattern does not get
 * 	a new line, it has to be parsed, when a pattern starts and ends
 *
 * @param pattern_split		The string containing a few pattern
 *
 * @param tokens			The allowed tokens which can be used to seperate patterns in a line
 */
std::vector<std::string> patternset::SplitString(std::string pattern_split, char* tokens) {
	std::vector<std::string> patternset;
	std::string tmp = "";
	bool flag_token = false;

	for (unsigned int i = 0; i < pattern_split.length(); i++) {
		for (unsigned int j = 0; j < strlen(tokens); j++) {
			if (pattern_split[i] == tokens[j]) {
				flag_token = true;
			}
		}
		if (flag_token) {								/*token found, which means in one line more patterns  --> start new pattern, save last pattern*/
			patternset.push_back(tmp);
			tmp = "";
		}
		else {
			tmp = tmp + pattern_split[i];					/*concatenating patternparts*/
		}
		flag_token = false;
	}
	patternset.push_back(tmp);
	return patternset;
}

/**
 * Validates a pattern, if it contains only pattern symbols and seperation tokens
 *
 * @param pattern_form		The pattern which has to be investigate for symbols and tokens
 *
 * @return 			Returns true if this one pattern is in right format, false else
 */
bool patternset::ValidatePatternsFormat(std::string pattern_form) {
	bool flag = true;

	for (unsigned int i = 0; i < pattern_form.length(); i++) {			/*allowed tokens in patternformat, also separating tokens*/
		if (pattern_form[i] != '1' && pattern_form[i] != '0' && pattern_form[i] != ' ' && pattern_form[i] != ',' && pattern_form[i] != '.' && pattern_form[i] != ';') {
			flag = false;
			SecureMessage("format", -1);
			break;
		}
	}
	if (pattern_form[0] != '1' || pattern_form[pattern_form.length() - 1] != '1') {
		flag = false;
		SecureMessage("startend", -1);
	}

	return flag;
}

/**
 * Validates a pattern set, if all patterns have the same length and weight
 *
 * @return 			Returns true if this one pattern is in right format, false else
 */
bool patternset::ValidatePatternConditions() {
	int com_weight, com_size;
	bool condition;

	com_size = (int) string_pat.size();
	com_weight = PatternWeight(string_pat[0]);					/*as fixpoint saving first patternweight*/
	condition = true;

	for (int i = 1; i < com_size; i++) {
		if (PatternWeight(string_pat[i]) != com_weight) {
			SecureMessage("weight", i);
			condition = false;
		}
	}
	return condition;
}

/**
 * Estimates the weight of a pattern
 *
 * @param pattern_str		The pattern which has to be investigate for the weight
 *
 * @return 			Returns true if this one pattern is in right format, false else
 */
int patternset::PatternWeight(std::string pattern_wght) {
	int str_weight;
	int str_length;

	str_weight = 0;
	str_length = (int) pattern_wght.length();

	for (int i = 0; i < str_length; i++) {
		if (pattern_wght[i] == '1') {
			str_weight++;
		}
	}
	return str_weight;
}

/**
 * Estimates the weight of a pattern
 *
 * @param pattern_str		The pattern which has to be investigate for the weight
 *
 * @return 			Returns true if this one pattern is in right format, false else
 */
int* patternset::PatternLength(std::vector<std::string> pattern_length) {
	int* pat_length;
	int minp, maxp;

	pat_length = new int[2];
	minp = (int) pattern_length[0].length();
	maxp = (int) pattern_length[0].length();

	for (unsigned int i = 1; i < pattern_length.size(); i++) {
		if ((int)pattern_length[i].length() > maxp) {
			maxp = (int) pattern_length[i].length();
		}
		if ((int)pattern_length[i].length() < minp) {
			minp = (int) pattern_length[i].length();
		}
	}
	pat_length[0] = minp;
	pat_length[1] = maxp;
	return pat_length;
}

/**
 * This function will create for our lenghts intervall
 * for each pattern a specific length!.
 *
 */
void patternset::CreateLengths(){
	double step, leng2;
	int leng, one;

	one = 0;
	leng = length[1];
	if(length[0] != length[1]){
		leng--;
	}
	if(length[0] == weight){
		one = 1;
	}
	std::uniform_int_distribution<int> pat_leng(length[0]+one, leng);

	if(randpatleng){
		for (int i = 0; i < size-1; i++) {
			if (i < (size / 4)) {
				lengths.push_back(length[1]);
			}
			else {
				lengths.push_back(pat_leng(generator));
			}
		}
		lengths.push_back(length[0]);
	}
	else{
		step = ((double)(length[1]-length[0]+1)/(size-1));
		leng = length[0]+one;
		leng2 = (double) leng;
		if(step == 0){
			step++;		
		}
		for(int i = one; i < size; i++){
			if(leng > length[1]){
				leng = length[1];
			}
			lengths.push_back(leng);
			leng2 += step;
			leng = (int)leng2;
		}
		if(one==1){				/*if one pattern has to be without don't care positions*/
			lengths.push_back(length[0]);
		}
	}
	std::sort(lengths.begin(),lengths.end());
}

/**
 * Creates random a set of pattern. For convention a pattern has to start and end with  '1'
 * 	Reason 10010 ~ 1001
 *
 * @return 			Returns a randomly created pattern set
 */
void patternset::CreateRandomPattern() {
	std::vector<uint64_t> tmp;
	uint64_t prototype, position, one;
	int counter;
	
	bool flag;

	one = 1;

	std::sort(lengths.rbegin(),lengths.rend());

	for (int i = 0; i < size; i++) {
		flag = true;
		while(flag){
			flag = false;
			prototype = 1;
			prototype = ((prototype << (lengths[i] - 1))| 1);					/*Prototype, first and last position is a match!*/
			std::uniform_int_distribution<int> distribution(1, lengths[i] - 2);			/*better random generater than time, uses likelihood for evenly random distribution*/
			counter = 2;
			while (counter < weight) {
				position = (uint64_t) distribution(generator);
				if ((prototype & (one << position)) == 0) {					/*have fun and fill with '1' randomly :D */
					counter++;
					prototype |= (one << position);
				}
			}
			for(int j = 0; j < i; j++){
				if(prototype == tmp[j]){
					flag = true;
				}
			}
			if(!flag){
				tmp.push_back(prototype);
			}
		}
	}
	for(int i = 0; i < size; i++){
		pattern_set.push_back(tmp[i]);
		string_pat.push_back(BitString(tmp[i]));
	}
	tmp.clear();
}


/**
 * Scans if there is another pattern in the same format
 *
 * @return returns boolean if there is another same pattern
 */
bool patternset::UniqPattern(int number) {
	bool uniq = true;
	for (int i = 0; i < size; i++) {
		if (number != i) {
			if (pattern_set[i] == pattern_set[number]) {
				uniq = false;
			}
		}
	}
	return uniq;
}

/*--Func-Change--------------------------------------------------------------*/
/**
 * Changes two different positions ('1' and '0') in a specific bit-pattern
 * Start and end are excluded
 *
 * @param number
 * 		The pattern which has to be modified
 */
void patternset::ChangeBits(int number) {
	int pos1, pos0;
	bool flag = true;

	if(lengths[number] == weight){		/*only possible, if lenghts[size-1] == weight*/
		number--;
	}

	while (flag && improve) {
		flag = false;
		pos1 = SymbolRandPos(number, '1');
		pos0 = SymbolRandPos(number, '0');
		pattern_set[number] = ChangeBitPos(number, pos1, pos0);
		if (!UniqPattern(number)) {
			flag = true;
		}
	}
}

/**
 * Exchange of two positions; changing match to don't care and vice versa.
 *
 * @param pos			The pattern index of the modifying pattern
 *
 * @param pos_one		Matchposition which has to be a don't care
 *
 * @param pos_zero		Don'tCare position which has to be a match
 * 
 * @return			The modified pattern as an 64-Bit integer
 */
uint64_t patternset::ChangeBitPos(int pos, int pos_one, int pos_zero) {
	uint64_t diff, one_mask, mask, p;

	p = pattern_set[pos];
	diff = 1;
	diff = diff << pos_one;
	one_mask = 1;
	
	one_mask = one_mask << ((int)log2(p)+1);
	one_mask--;

	mask = one_mask ^ diff;
	p = p & mask;

	diff = 1;
	diff = diff << pos_zero;
	p = p | diff;

	return p;
}

/**
 * Returns for a specific symbol ('1', '0') a random chosen position
 * 
 * @param number		The pattern index of the modifying pattern
 *
 * @param symb			The specific symbol, either '1' or '0'
 *
 * @return			A random chosen position in the pattern of the specific symbol
 */
int patternset::SymbolRandPos(int number, char symb) {
	std::vector<int> positions;
	int pos;

	positions = GetSymbol(number, symb);

	std::uniform_int_distribution<int> distribution(0, (int) positions.size() - 1);

	pos = distribution(generator);
	pos = positions[pos];

	positions.clear();
	return pos;
}

/**
 * Investigates for a specific symbol all positions of a pattern.
 * For a match, the first and last positions are excluded!
 * 
 * @param number		The pattern index of the modifying pattern
 *
 * @param symb			The specific symbol, either '1' or '0'
 *
 * @return 			All positions of the pattern of the specific symbol
 */
std::vector<int> patternset::GetSymbol(int number, char symb){
	std::vector<int> positions;
	int lgth;
	uint64_t tmp, chr;

	tmp = pattern_set[number];
	lgth = (int)(log2(tmp));

	chr = 0;
	if (symb == '1') {
		chr = 1;
	}
 
	for (int i = 1; i < lgth; i++) {
		if (((tmp >> i) & 1) == chr) {
			positions.push_back(i);
		}
	}
	return positions;		
}


/*---stuff-------------------------------------------------------------------*/
/**
* Method to print the current pattern
*/
void patternset::Print() {
	ToString();
	for (int i = 0; i < size; i++) {
		std::cout << string_pat[i] << std::endl;;
	}
}

/**
 * Saves the current 64-Bit integer in the string pat as printable strings
 */
void patternset::ToString() {
	for (int i = 0; i < size; i++) {
		string_pat[i] = BitString(pattern_set[i]);
	}
}

/**
 * Turns a string into a 64-Bit Integer as binary pattern 
 * 
 * @param p 			The number of the pattern in the string pattern set
 *
 * @return 			The binary representation of the input string
 */
uint64_t  patternset::ToBit(std::string p) {
	std::string tmp;
	uint64_t val, d;
	val = 0;
	d = 1;

	for (unsigned int i = 0; i < p.length(); i++) {
		if (p[i] == '1') {
			val += d << (p.length() - (i + 1));
		}
	}
	return val;
}

/**
 * Changes a pattern in 64-Bit integer format into a printable string
 *
 * @param p			The pattern in 64-bit integer format
 *
 * @return 			The printable string of the input pattern
 */
std::string patternset::BitString(uint64_t p) {
	std::string bit_string;
	int lgth;
	lgth = (int) log2(p);

	bit_string = "";

	for (int i = lgth; i >= 0; i--) {
		if (((p >> i) & 1) == 1) {
			bit_string += "1";
		}
		else {
			bit_string += "0";
		}
	}
	bit_string += "\0";
	return bit_string;
}


/**
 * Determines the maximum number of patterns with weight and length
 *
 *@param weight			Complete pattern weight-2, start and end have to be match and do not change
 *
 *@param length			Complete pattern lengt-2, start and end have to be match and do not change
 *
 *@return 			max number of possible pattern
 */
double patternset::MaxNumberPattern(int p_weight, int p_length) {
	double tmpa = Faculty(p_length);
	double tmpb = Faculty(p_length - p_weight);
	double tmpc = Faculty(p_weight);
	double tmp = tmpa / (tmpb*tmpc);
	return tmp;
}

/**
 * Calculates the faculty
 *
 *@param value			Calculating the faculty for value
 *
 *@return 			faculty(value) = 'value!'
 */
double patternset::Faculty(int value) {
	double tmp;
	tmp = 1;
	for (int i = 1; i <= value; i++) {
		tmp *= i;
	}
	return tmp;
}

/**
 * A Method to collect all errormessages. Just easier for programmer to change
 *  	the text or extend.
 *
 * @param errmsg			Due to a few possible errormessages, this is the option, which has to be printed.
 *
 * @param pos			The position of the incorrect patterns.
 */
void patternset::SecureMessage(std::string errmsg, int pos) {
	/*if (errmsg == "wrongindex") {
		std::cerr << "ERROR! Pattern " << pos << " does not exist... do nothing\n" << std::endl;
		return;
	}
	if (errmsg == "file") {
		std::cerr << "ERROR! Pattern file \'" << pattern_file << "\' could not be found!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if (errmsg == "empty") {
		std::cerr << "ERROR! File \'" << pattern_file << "\' is an empty file!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if (errmsg == "pattern") {
		std::cerr << "ERROR! Patternconditions from pattern file were not correct (different weight or length)!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if (errmsg == "nkl") {
		std::cerr << "ERROR! Wrong values for weight, pattern number or pattern length!" << std::endl;
		std::cerr << "Return to default values\n" << std::endl;
		return;
	}
	if (errmsg == "weight_pat") {
		std::cerr << "ERROR! Weight of a pattern cannot under 2!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if (errmsg == "format") {
		std::cerr << "FORMAT-ERROR: Pattern containes illegal characters!" << std::endl;
		std::cerr << "Allowed characters: '0','1' for pattern; ','|'.'|';'|' ' to seperate patterns" << std::endl;
		std::cerr << "Go on to next pattern.\n" << std::endl;
		return;
	}
	if (errmsg == "max_number_pattern") {
		std::cerr << "The number of patterns is too high for your configuration and will be reset!" << std::endl;
		return;
	}
	if (errmsg == "size") {
		std::cerr << "ERROR! The number of patterns has to be greater than 0!\n" << std::endl;
		return;
	}
	if (errmsg == "weight") {
		std::cerr << "ERROR! By comparing with the first pattern, the " << pos + 1 << ". pattern has a different weight!\n" << std::endl;
		return;
	}
	if(errmsg == "length"){
		std::cerr << "ERROR! The patternlengths has to be greater or equal the weight, greater than 2 and below 63!" << std::endl;
		return;
	}
	if (errmsg == "startend") {
		std::cerr << "FORMAT-ERROR: Pattern has to start and end with a match position '1' !\n" << std::endl;
		return;
	}
	if (errmsg == "noimprove") {
		std::cerr << "Using your pattern conditions it is not sensible to improve your pattern, sorry!" << std::endl;
		std::cerr << "Deactivating improve mode\n" << std::endl;
		return;
	}
	if (errmsg == "update"){
		std::cerr << "\n####IMPROTANT####\nDue to some configuration errors, your submitted parameters have been updated!" << std::endl;
		return;
	}*/
}


/*---Set-&-GetFunc-----------------------------------------------------------*/
/**
 * Changes to silent output, nothing will be printed automatically!
 * Errors are going to be printed!
 */
void patternset::Silent(){
	silent = true;
}

/**
 * Changes the patternlengths from evenly distributed to random pattern lengths.
 * One quarter has maximum length, the rest randomly distributed between 
 *	maximum and minimum length.
 */
void patternset::RandPatLength(){
	if(!randpatleng){
		randpatleng = true;
		ReInitPattern();
	}
}

/**
 * Returns the weight of each pattern, the match positions
 *
 * @return 		returns weight
 */
int patternset::GetWeight() {
	return weight;
}

/**
 * Returns the amount of patters; number of patterns
 *
 * @return 		returns size
 */
int patternset::GetSize() {
	return size;
}

/**
 * Returns the length of each pattern
 *
 * @return 		returns length
 */
int* patternset::GetLength() {
	return length;
}

/**
 * Returns the mean of all pattern lengths
 * 
 * @return 		The mean of all lengths;
 */
int patternset::GetLengthMean(){
	return length_mean;
}

/**
 * Returns a boolean, if an improvement is possible
 * Imagine, a improvment for pattern '101' is not possible!
 * 
 * @return 		The boolean, if this improvement is possible
 */
bool patternset::GetImprove() {
	return improve;
}

/**
 * Returns a boolean, if due to some configuration errors
 *	values set by the user had to be updated to better configuration.
 * Imagine, the weight was above the pattern length.
 *
 * @return 		The boolean, if something had been updated.
 */
bool patternset::GetUpdate() {
	return update;
}


/**
 * Returns a specific pattern of the pattern set in 64-Bit integer format
 *
 * @return 		specific string in 64-Bit integer format
 */
uint64_t patternset::GetPattern(int number) {
	if (number >= size) {
		SecureMessage("wrongindex", number);
		return 0;
	}
	else {
		return pattern_set[number];
	}
}

/**
 * Returns complete pattern set
 *
 * @return 		complete pattern set
 */
std::vector<std::string> patternset::GetStringPattern() {
	ToString();
	return string_pat;
}

/**
 * Returns only one pattern as a string from the current patternset.
 *
 * @param number	The index of the chosen pattern
 *
 * @return 		The pattern as a printable string
 */
std::string patternset::GetString(int number){
	if (number >= size) {
		SecureMessage("wrongindex", number);
		return NULL;
	}
	else {
		ToString();
		return string_pat[number];
	}
}

/**
 * Changes a pattern in 64-Bit integer format
 *
 * @param number	The index of the pattern which will be changed
 *
 * @param patt		The new pattern
 */
void patternset::SetPattern(int number, uint64_t patt){
	if (number >= size) {
		SecureMessage("wrongindex", number);
	}
	else {
		pattern_set[number] = patt;
	}
}
