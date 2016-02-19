/** 
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * extended key object file
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
#include "extkey.h"

/*---Variables---------------------------------------------------------------*/


/*===Main-Part===============================================================*/
/*---Constructor-& Init------------------------------------------------------*/
extkey::extkey(){
	extkey(1,1);
}

extkey::extkey(int pos, double value){
	this->pos = pos;
	this->value = value;
}


/*---Functions---------------------------------------------------------------*/
int extkey::GetPos(){
	return pos;
}

double extkey::GetValue(){
	return value;
}

void extkey::SetValue(double value){
	this->value = value;
}

bool extkey::operator < (extkey const& p) const{
		return (value < p.value);
}

bool extkey::operator <= (extkey const& p) const{
		return (value <= p.value);
}

bool extkey::operator > (extkey const& p) const{
		return (value > p.value);
}

bool extkey::operator >= (extkey const& p) const{
		return (value >= p.value);
}

bool extkey::operator == (extkey const& p) const{
		return (value == p.value);
}

bool extkey::operator || (extkey const& p) const{
		return (value || p.value);
}

bool extkey::operator != (extkey const& p) const{
		return (value != p.value);
}
