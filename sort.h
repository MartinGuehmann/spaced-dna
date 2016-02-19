/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#ifndef SORT_H_
#define SORT_H_

#include <stdio.h>
 #include <string.h>
#include <getopt.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <ctype.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <tr1/unordered_map>
#include <map>
#include <limits>
#include <unordered_map>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iomanip>
#include <random>
#include <bitset>

#define DNA_BUCKETS 256
#define PROT_BUCKETS 26
#define EU 0
#define JS 1
#define EV 2

using namespace std;

template<typename uint>
struct word{ 
	unsigned long long key;
	unsigned short seq;
	bool operator<( const word& val ) const { 
		return key < val.key; 
	}
};

template<typename uint>
struct sequence{ 
	uint headerStart;
	uint headerEnd;
	uint seqStart;
	uint seqEnd;
	uint seqWordEnd;
	uint length=0;
	uint revCompStart;
	uint revCompEnd;
	uint revCompWordEnd;
	vector<double> frequencies={0,0,0,0};
	vector<double> frequenciesRevComp={0,0,0,0};
};

void getFirstBits(vector<unsigned char>&pattern, unsigned char* firstBits, unsigned long long pos, unsigned long length, unsigned char* seqData, int weight){
	unsigned char alphabet[256]={};
	alphabet['A']=0;alphabet['C']=1;alphabet['G']=2;alphabet['T']=3; // mapping A -> 00, C -> 01 etc. without overhead
	unsigned char w;
	// mapping like:
	// ACAG -> w = 00 01 00 10 
	// uses one byte for first bits, thus max weight = 4
	unsigned char p0,p1,p2,p3;
	p0=pattern[0];p1=pattern[1];p2=pattern[2];p3=pattern[3];
	if(weight>=4){
		for(int i=0; i<length;i++){
			w=0;
			unsigned char* tmp=&seqData[i+pos];
			w|=alphabet[*(tmp+p0)] << 6;
			w|=alphabet[*(tmp+p1)] << 4;
			w|=alphabet[*(tmp+p2)] << 2;
			w|=alphabet[*(tmp+p3)];
			firstBits[i]=w;
		}
	}
	else if (weight==3){
		unsigned char p0,p1,p2;
		p0=pattern[0];p1=pattern[1];p2=pattern[2];
		for(int i=0; i<length;i++){
			w=0;
			unsigned char* tmp=&seqData[i+pos];
			w|=alphabet[*(tmp+p0)] << 4;
			w|=alphabet[*(tmp+p1)] << 2;
			w|=alphabet[*(tmp+p2)];
			firstBits[i]=w;
		}
	}
	else if(weight==2){
		unsigned char p0,p1;
		p0=pattern[0];p1=pattern[1];
		for(int i=0; i<length;i++){
			w=0;
			unsigned char* tmp=&seqData[i+pos];
			w|=alphabet[*(tmp+p0)] << 2;
			w|=alphabet[*(tmp+p1)];
			firstBits[i]=w;
		}
	
	}
	else{
		unsigned char p0;
		p0=pattern[0];
		for(int i=0; i<length;i++){
			unsigned char* tmp=&seqData[i+pos];
			w=alphabet[*(tmp+p0)];
			firstBits[i]=w;
		}
	}
}

template<typename uint>
void writeDmat(vector<vector<double> > dmat, vector<sequence<uint> >& sequences, string filename, unsigned char* seqData){
	ofstream outfile;
	outfile.open(filename.c_str());
	outfile << sequences.size() << endl;
	for (int i = 0; i < sequences.size(); i++) {
		for(uint k=0; k<10;k++){
			if( k+sequences[i].headerStart <= sequences[i].headerEnd )
				outfile << seqData[k+sequences[i].headerStart];
			else
				outfile << " ";
		}
		outfile << " ";
     	for (int j = 0; j < sequences.size(); j++) {
			if (i > j) 
	    			outfile << setprecision(12) << dmat[i][j] << "  ";
			else if(j>i)
				outfile << setprecision(12) << dmat[j][i] << "  ";
			else
					outfile << setprecision(12) << "0" << "  ";
     	}
      		outfile << endl;
	}
	//outfile << dmat[1][0] << endl;
}

template<typename uint>
void spacedDNA(vector<string>& patternSet, unsigned char* seqData, uint n, int distance, int threads, int weight, int dontCare, bool revComp, string output){
	int ell=dontCare+weight-1;
	#ifdef _OPENMP
	omp_set_dynamic(0);
	omp_set_num_threads(threads);
	#endif
	unsigned char alphabet[256]={};
	alphabet['A']=0;alphabet['C']=1;alphabet['G']=2;alphabet['T']=3;
	vector<sequence<uint> > sequences; 
	int seqNum;
	uint totalSeqLength=0;
	for(uint i=0,j=0; i<n;i++){
		if(seqData[i]==';'){ //remove comments
			while(seqData[++i]!='\n');
		}else if(seqData[i]=='>'){//read header
			sequence<uint> tmpSeq={};
			seqData[j++]=seqData[i];
			tmpSeq.headerStart=j;
			totalSeqLength++;
			if(sequences.size()>0){
				sequences.back().seqEnd=j-1;
				sequences.back().seqWordEnd=j-(weight+dontCare);
			}
			while(seqData[++i]!='\n'){
				seqData[j++]=seqData[i];
				totalSeqLength++;
			}
			tmpSeq.headerEnd=j-1;
			tmpSeq.seqStart=j;
			sequences.push_back(tmpSeq);
		}else if(isalpha(seqData[i])){
			char c=toupper(seqData[i]);
			switch(c){
				case 'A': 
					seqData[j++]='A';
					sequences.back().frequencies[0]++;
					sequences.back().length++;
					break;
				case 'C': 
					seqData[j++]='C';
					sequences.back().frequencies[1]++;
					sequences.back().length++;
					break;
				case 'G': 
					seqData[j++]='G';
					sequences.back().frequencies[2]++;
					sequences.back().length++;
					break;
				case 'T': 
					seqData[j++]='T';
					sequences.back().frequencies[3]++;
					sequences.back().length++;
					break;
				default: 
					seqData[j++]='N';
					sequences.back().length++;
					break;
			}
			totalSeqLength++;
		}
	}
	sequences.back().seqEnd=totalSeqLength;
	sequences.back().seqWordEnd=totalSeqLength-(weight+dontCare)+1;
	seqNum=sequences.size();
	uint cnt=0;
	for(uint i=0; i<seqNum;i++){
		uint start=sequences[i].seqStart;
		uint end=sequences[i].seqEnd-1;
		sequences[i].revCompStart=cnt+totalSeqLength;
		sequences[i].revCompEnd=cnt+totalSeqLength+sequences[i].length;
		sequences[i].revCompWordEnd=sequences[i].revCompEnd-(weight+dontCare)+1;
		sequences[i].frequenciesRevComp[0]=sequences[i].frequencies[3];
		sequences[i].frequenciesRevComp[1]=sequences[i].frequencies[2];
		sequences[i].frequenciesRevComp[2]=sequences[i].frequencies[1];
		sequences[i].frequenciesRevComp[3]=sequences[i].frequencies[0];
		for(uint j=end; j>=start ;j--){
			switch(seqData[j]){
				case 'A': seqData[cnt++ + totalSeqLength]='T';break;
				case 'C': seqData[cnt++ + totalSeqLength]='G';break;
				case 'G': seqData[cnt++ + totalSeqLength]='C';break;
				case 'T': seqData[cnt++ + totalSeqLength]='A';break;
				case 'N': seqData[cnt++ + totalSeqLength]='N';break;
				default: exit(-1);
			}
		}
	}
	for(int i=0; i< sequences.size();i++){
		for(int j=0; j<4;j++){
			sequences[i].frequencies[j]=sequences[i].frequencies[j]/sequences[i].length;
			sequences[i].frequenciesRevComp[j]=sequences[i].frequenciesRevComp[j]/sequences[i].length;
		}
	}
	totalSeqLength*=2;
	vector<unsigned char> firstBits(totalSeqLength);
	vector<uint> prefixSum(256,0);
	vector<unsigned char> matchPos;
	vector<unsigned char> dCarePos;
	vector<vector <vector <double> > > dmat(threads, vector< vector<double> >(seqNum, vector<double>(seqNum,0)));
	vector<vector <double> > dmatFinal(seqNum, vector<double>(seqNum,0));
	vector<vector<double> > row(threads, vector<double>(seqNum*2));
	for(int p = 0; p < patternSet.size(); p++ ){
		for(int k = 0; k < patternSet[p].length(); k++){
			if (patternSet[p][k] == '1')
				matchPos.push_back(k);
			else
				dCarePos.push_back(k);
		}
		getFirstBits(matchPos, &firstBits[0], 0, totalSeqLength, seqData, weight);
		fill(prefixSum.begin(), prefixSum.end(), 0);

		for(int i=0;i<seqNum;i++){
			for(uint j=sequences[i].seqStart; j<sequences[i].seqWordEnd;j++){
				++prefixSum[(firstBits[j])];	
			}
			for(uint j=sequences[i].revCompStart; j<sequences[i].revCompWordEnd;j++){
				++prefixSum[(firstBits[j])];	
			}
		}
		#pragma omp parallel for schedule(runtime)
		for(unsigned int k=0; k<DNA_BUCKETS; k++){
			vector<word<uint> > words(prefixSum[k]);
			double cnt=0;	
			unsigned long long w;
			unsigned char bits;
			for(int i=0;i<seqNum;i++){
				for(uint j=sequences[i].seqStart; j<sequences[i].seqWordEnd;j++){
					if(firstBits[j]==k){
						w=0;
						bits=weight*2-2;
						unsigned char c;
						bool correctWord=true;
						for(int o=0; o<matchPos.size();o++){
							uint pos=matchPos[o];
							c= seqData[j+pos];
							if(c=='N')
								correctWord=false;
							w|=(unsigned long long) alphabet[c] << bits;
							bits-=2;
						}
						if(correctWord){
							words[cnt].key=w;
							words[cnt++].seq=i;
						}
					}
				}
				if(revComp){
					for(uint j=sequences[i].revCompStart; j<sequences[i].revCompWordEnd;j++){
						if(firstBits[j]==k){
							w=0;
							bits=weight*2-2;
							unsigned char c;
							bool correctWord=true;
							for(int o=0; o<matchPos.size();o++){
								uint pos=matchPos[o];
								c= seqData[j+pos];
								if(c=='N')
									correctWord=false;
								w|=(unsigned long long) alphabet[c] << bits;
								bits-=2;
							}
							if(correctWord){
								words[cnt].key=w;
								words[cnt++].seq=i+seqNum;
							}
						}		
					}	
				}		
			}
			std::sort(words.begin(),words.begin()+cnt);
			int thread=0;
			#ifdef _OPENMP
			thread=omp_get_thread_num();
			#endif
			if(distance==EU){
				fill(row[thread].begin(), row[thread].end(), 0);
				for(uint i=0; i<cnt;){
					unsigned long long key_t=words[i].key;
					row[thread][words[i].seq]++;
					while(++i<cnt &&key_t==words[i].key){
						row[thread][words[i].seq]++;
					}
						for(int i=0;i<seqNum;i++){
							for(int j=0; j<i; j++){
								double t1;
								if(revComp)
									t1=abs(row[thread][i]-(row[thread][j]+row[thread][j+seqNum]));
								else
									t1=abs(row[thread][i]-row[thread][j]);
								dmat[thread][i][j]+=t1*t1;
							}	
						}
					fill(row[thread].begin(), row[thread].end(), 0);	
				}
			}
			else if(distance==EV){
				fill(row[thread].begin(), row[thread].end(), 0);
				for(uint i=0; i<cnt;){
					unsigned long long key_t=words[i].key;
					row[thread][words[i].seq]++;
					while(++i<cnt &&key_t==words[i].key){
						row[thread][words[i].seq]++;
					}
					for(int i=0;i<seqNum;i++){
						for(int j=0; j<i; j++){
							if(revComp)
								dmat[thread][i][j]+=min(row[thread][i],row[thread][j]+row[thread][j+seqNum]);
							else
								dmat[thread][i][j]+=min(row[thread][i],row[thread][j]);
						}	
					}
					fill(row[thread].begin(), row[thread].end(), 0);	
				}
			}
			else if(distance==JS){
				fill(row[thread].begin(), row[thread].end(), 0);
				for(uint i=0; i<cnt;){
					unsigned long long key_t=words[i].key;
					row[thread][words[i].seq]++;
					while(++i<cnt &&key_t==words[i].key){
						row[thread][words[i].seq]++;
					}
					for(int i=0; i<seqNum; i++){
						if(revComp){
							row[thread][i+seqNum]=(row[thread][i]+row[thread][i+seqNum])/((sequences[i].seqWordEnd-sequences[i].seqStart)*2);
							row[thread][i]/=sequences[i].seqWordEnd-sequences[i].seqStart;
						}
						else
							row[thread][i]/=sequences[i].seqWordEnd-sequences[i].seqStart;
					}
					for(int i=0;i<seqNum;i++){
						double row_i=row[thread][i];
						for(int j=0; j<i; j++){
							double m;
							double row_j;
							if(revComp)
								row_j=row[thread][j+seqNum];
							else
								row_j=row[thread][j];
							m=0.5*(row_i+row_j);
							if(m==0)
								continue;
							if(row_i!=0 &&row_j!=0){
								dmat[thread][i][j]+= log2(row_i/m ) * row_i;
								dmat[thread][j][i]+= log2(row_j/m) * row_j;
							}
							else if(row_i!=0)
								dmat[thread][i][j]+= log2(row_i/m) * row_i;
							else if(row_j!=0)
								dmat[thread][j][i]+= log2(row_j/m) * row_j;
						}	
					}
					fill(row[thread].begin(), row[thread].end(), 0);	
				}
			}
		}
		matchPos.clear();
		dCarePos.clear();
		if(distance==EU || distance==JS){
			for(int p=1; p<threads; p++){
				for(int i=0;i<seqNum;i++){
					for(int j=0; j<seqNum; j++){
						dmat[0][i][j]+=dmat[p][i][j];
						dmat[p][i][j]=0;
					}
				}
			}
			if(distance==EU){
				#ifdef _OPENMP
				#pragma omp parallel for schedule(static)
				#endif
				for(int i=0;i<seqNum;i++){
					for(int j=0; j<i; j++)
						dmat[0][i][j]=sqrt(dmat[0][i][j]);
				}
				for(int i=0;i<seqNum;i++){
					for(int j=0; j<i; j++){
						dmatFinal[i][j]+=dmat[0][i][j];
						dmat[0][i][j]=0;
					}
				}
			}
			if(distance==JS){
				#ifdef _OPENMP
				#pragma omp parallel for schedule(static)
				#endif
				for(int i=0;i<seqNum;i++){
					for(int j=0; j<i; j++){
						dmat[0][i][j]=(0.5*dmat[0][i][j]+0.5*dmat[0][j][i]);
					}
				}
				for(int i=0;i<seqNum;i++){
					for(int j=0; j<i; j++){
						dmatFinal[i][j]+=dmat[0][i][j];
						dmat[0][i][j]=0;
						dmat[0][j][i]=0;
					}
				}
			}
		}
	}
	if(patternSet.size()>1&&distance!=EV){
		for(int i=0;i<seqNum;i++){
			for(int j=0; j<i; j++)
				dmatFinal[i][j]/=patternSet.size();
		}
	}
	else if(distance==EV){
		for(int p=1; p<threads; p++){
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<seqNum; j++){
					dmat[0][i][j]+=dmat[p][i][j];
				}
			}
		}
		#ifdef _OPENMP
  		#pragma omp parallel for schedule(static)
  		#endif
		for(int i=0;i<seqNum;i++){
			for(int j=0; j<i; j++){
				double min=std::min(sequences[i].length,sequences[j].length)-ell;
				double max=std::max(sequences[i].length,sequences[j].length)-ell;
				double q=0;
				for(int k=0; k<4;k++){
					if(revComp)
 						q+=((sequences[i].frequencies[k]+sequences[i].frequenciesRevComp[k])*0.5)*((sequences[j].frequencies[k]+sequences[j].frequenciesRevComp[k])*0.5);
					else
						q+=sequences[i].frequencies[k]*sequences[j].frequencies[k];
				}	
				double valueUnderRoot;
				if(revComp)
					valueUnderRoot=dmat[0][i][j]/(patternSet.size()*min)-2*max*pow(q,weight);
				else
					valueUnderRoot=dmat[0][i][j]/(patternSet.size()*min)-max*pow(q,weight);
				if(valueUnderRoot>=0){
					double p = pow(valueUnderRoot, (1.0 / weight));
					dmat[0][i][j]=-0.75 * log((4.0/3.0)*p-(1.0/3.0));
				}
				else{
					dmat[0][i][j]=1.2;
				}
			}
		}
		for(int i=0;i<seqNum;i++){
			for(int j=0; j<i; j++)
				dmatFinal[i][j]+=dmat[0][i][j];
		}
	}
	writeDmat<uint>(dmatFinal,sequences, output, seqData);
}

template<typename uint>
void spacedProt(vector<string>& patternSet, unsigned char* seqData, uint n, int distance, int threads, int weight, int dontCare, string output){
	#ifdef _OPENMP
	omp_set_dynamic(0);
	int ell=dontCare+weight-1;
	omp_set_num_threads(threads);
	#endif
	unsigned char alphabet[256]={};
	alphabet['A']=0;alphabet['B']=1;alphabet['C']=2;alphabet['D']=3;alphabet['E']=4;alphabet['F']=5;alphabet['G']=6;alphabet['H']=7;alphabet['I']=8;alphabet['J']=9;
	alphabet['K']=10;alphabet['L']=11;alphabet['M']=12;alphabet['N']=13;alphabet['O']=14;alphabet['P']=15;alphabet['Q']=16;alphabet['R']=17;alphabet['S']=18;alphabet['T']=19;
	alphabet['U']=20;alphabet['V']=21;alphabet['W']=22;alphabet['X']=23;alphabet['Y']=24;alphabet['Z']=25;
	vector<sequence<uint> > sequences; 
	int seqNum;
	uint totalSeqLength=0;
	for(uint i=0,j=0; i<n;i++){
		if(seqData[i]==';'){ //remove comments
			while(seqData[++i]!='\n');
		}else if(seqData[i]=='>'){//read header
			sequence<uint> tmpSeq={};
			seqData[j++]=seqData[i];
			tmpSeq.headerStart=j;
			totalSeqLength++;
			if(sequences.size()>0){
				sequences.back().seqEnd=j-1;
				sequences.back().seqWordEnd=j-(weight+dontCare);
			}
			while(seqData[++i]!='\n'){
				seqData[j++]=seqData[i];
				totalSeqLength++;
			}
			tmpSeq.headerEnd=j-1;
			tmpSeq.seqStart=j;
			sequences.push_back(tmpSeq);
		}else if(isalpha(seqData[i])){
			char c=toupper(seqData[i]);
			seqData[j++]=c;
			sequences.back().length++;
			totalSeqLength++;
		}
	}
	sequences.back().seqEnd=totalSeqLength;
	sequences.back().seqWordEnd=totalSeqLength-(weight+dontCare)+1;
	seqNum=sequences.size();
	vector<uint> prefixSum(PROT_BUCKETS,0);
	vector<unsigned char> matchPos;
	vector<unsigned char> dCarePos;
	vector<vector <vector <double> > > dmat(threads, vector< vector<double> >(seqNum, vector<double>(seqNum)));
	vector<vector <double> > dmatFinal(seqNum, vector<double>(seqNum,0));
	vector<vector<double> > row(threads, vector<double>(seqNum));
	for(int p = 0; p < patternSet.size(); p++ ){
		for(int k = 0; k < patternSet[p].length(); k++){
			if (patternSet[p][k] == '1')
				matchPos.push_back(k);
			else
				dCarePos.push_back(k);
		}
		fill(prefixSum.begin(), prefixSum.end(), 0);
		for(int i=0;i<seqNum;i++){
			for(uint j=sequences[i].seqStart; j<sequences[i].seqWordEnd;j++){
				++prefixSum[alphabet[seqData[j]]];	
			}
		}
		#pragma omp parallel for schedule(runtime)
		for(unsigned int k=0; k<PROT_BUCKETS; k++){
			vector<word<uint> > words(prefixSum[k]);
			double cnt=0;	
			unsigned long long w;
			unsigned char bits;
			for(int i=0;i<seqNum;i++){
				for(uint j=sequences[i].seqStart; j<sequences[i].seqWordEnd;j++){
					if(alphabet[seqData[j]]==k){
						w=0;
						bits=weight*5-5;
						unsigned char c;
						for(int o=0; o<matchPos.size();o++){
							uint pos=matchPos[o];
							c= seqData[j+pos];
							w|=(unsigned long long) alphabet[c] << bits;
							bits-=5;
						}
						words[cnt].key=w;
						words[cnt++].seq=i;
					}
				}
			}
			std::sort(words.begin(),words.begin()+cnt);
			int thread=0;
			#ifdef _OPENMP
			thread=omp_get_thread_num();
			#endif
			if(distance==EU){
				fill(row[thread].begin(), row[thread].end(), 0);
				for(uint i=0; i<cnt;){
					unsigned long long key_t=words[i].key;
					row[thread][words[i].seq]++;
					while(++i<cnt &&key_t==words[i].key){
						row[thread][words[i].seq]++;
					}
					for(int i=0;i<seqNum;i++){
						for(int j=0; j<i; j++){
							double t1;
							t1=abs(row[thread][i]-row[thread][j]);
							dmat[thread][i][j]+=t1*t1;
						}	
					}
					fill(row[thread].begin(), row[thread].end(), 0);	
				}
			}
			else if(distance==JS){
				fill(row[thread].begin(), row[thread].end(), 0);
				for(uint i=0; i<cnt;){
					unsigned long long key_t=words[i].key;
					row[thread][words[i].seq]++;
					while(++i<cnt &&key_t==words[i].key){
						row[thread][words[i].seq]++;
					}
					for(int i=0; i<seqNum; i++){
						row[thread][i]/=sequences[i].seqWordEnd-sequences[i].seqStart;
					}
					for(int i=0;i<seqNum;i++){
						double row_i=row[thread][i];
						for(int j=0; j<i; j++){
							double m;
							double row_j;
							row_j=row[thread][j];
							m=0.5*(row_i+row_j);
							if(m==0)
								continue;
							if(row_i!=0 &&row_j!=0){
								dmat[thread][i][j]+= log2(row_i/m ) * row_i;
								dmat[thread][j][i]+= log2(row_j/m) * row_j;
							}
							else if(row_i!=0)
								dmat[thread][i][j]+= log2(row_i/m) * row_i;
							else if(row_j!=0)
								dmat[thread][j][i]+= log2(row_j/m) * row_j;
						}	
					}
					fill(row[thread].begin(), row[thread].end(), 0);	
				}
			}
		}
		matchPos.clear();
		dCarePos.clear();
		for(int p=1; p<threads; p++){
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<seqNum; j++){
					dmat[0][i][j]+=dmat[p][i][j];
				}
			}
		}
		if(distance==EU){
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<i; j++)
					dmat[0][i][j]=sqrt(dmat[0][i][j]);
			}
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<i; j++)
					dmatFinal[i][j]+=dmat[0][i][j];
			}
		}
		if(distance==JS){
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<i; j++){
					dmat[0][i][j]=(0.5*dmat[0][i][j]+0.5*dmat[0][j][i]);
				}
			}
			for(int i=0;i<seqNum;i++){
				for(int j=0; j<i; j++)
					dmatFinal[i][j]+=dmat[0][i][j];
			}
		}
	}
	if(patternSet.size()>1&&distance!=EV){
		for(int i=0;i<seqNum;i++){
			for(int j=0; j<i; j++)
				dmatFinal[i][j]/=patternSet.size();
		}
	}
	writeDmat<uint>(dmatFinal,sequences, output, seqData);
}

#endif 

