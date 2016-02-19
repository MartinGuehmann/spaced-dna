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

#include "sort.h"
#include "variance.h"

bool revComp=true;
int number = 5, dontcare = 15, weight = 14, threads = 15, distanceType = EV;
string patternFile="", output="DMat";

void printHelp(){
 string help = "Alignment-free sequence comparison tool using spaced word frequencies"
    "\nUsage: ./spaced [options] <sequence file>"
    "\nOptions:"        
    "\n\t -h: print this help and exit"
    "\n\t -o <file>: output file name (default: DMat)"
    "\n\t -k <integer>: pattern weight (default 14)"
    "\n\t -l <integer>: pattern don't care positions (default 15)"
    "\n\t -n <integer>: number of patterns (default: 5)"  
    "\n\t -f <file>: use patterns in <file> instead of random patterns"  
    "\n\t -t <integer>: numer of threads (default: 25 threads)"
    "\n\t -r: don't consider the reverse complement"
    "\n\t -d EU | JS | EV: change distance type to Euclidean, Jensen-Shannon, evolutionary distance (default: EV)" 
    "\n";
	cout << help << endl;
}

void parsePattern(string filePath, vector<string>& patternSet){
    ifstream file(filePath.c_str());
    string str; 
    while (getline(file, str)){
            patternSet.push_back(str);
        }
    file.close();
}

void parseParameters(int argc, char** argv){
    // too few?
    if (argc < 2) {
        printHelp(); // print help
        exit(0); // and exit
    }
    // only parameter -h?
    if (argc == 2)
        if (strcmp("-h", argv[1]) == 0) {
            printHelp();
            exit(0);
    	}
    for (int i = 1; i < argc - 1; i++) { // first argument is run dir, last should be sequence file
        // help & exit
        if (strcmp("-h", argv[i]) == 0) {
            printHelp();
            exit(0);
        }
        else if (strcmp("-r", argv[i]) == 0) {
            revComp=false;
        }
        // patterns that should be used
        else if (strcmp("-f", argv[i]) == 0) {
            patternFile = argv[++i];
            if (i == argc - 1) {
                cerr << "Flag -f needed an argument!" << endl;
                exit(-1);
            }
        }
        // output matrix file
        else if (strcmp("-o", argv[i]) == 0) {
            output = argv[++i];
            if (i == argc - 1) {
                cerr << "Flag -o needed an argument!" << endl;
                exit(-1);
            }            
        }
        // auto-generate pattern weight
        else if (strcmp("-k", argv[i]) == 0) {
            weight = atoi(argv[++i]);
            if (i == argc - 1) {
                cerr << "Flag -k needed an argument!" << endl;
                exit(-1);
            }                
        }
        // auto-generate pattern length
        else if (strcmp("-l", argv[i]) == 0) {
            dontcare = atoi(argv[++i]);
            if (i == argc - 1) {
                cerr << "Flag -l needed an argument!" << endl;
                exit(-1);
            }
        }
        // number of combined matrices
        else if (strcmp("-n", argv[i]) == 0) {
            number = atoi(argv[++i]);
            if (i == argc - 1) {
                cerr << "Flag -n needed an argument!" << endl;
                exit(-1);
            }            
        }
        // number of threads
        else if (strcmp("-t", argv[i]) == 0) {
            threads = atoi(argv[++i]);
            if (i == argc - 1) {
                cerr << "Flag -t needed an argument!" << endl;
                exit(-1);
            }            
        }
        // distance type
        else if (strcmp("-d", argv[i]) == 0) {
            i++;
            if (i == argc - 1) {
                cerr << "Flag -d needed an argument!" << endl;
                exit(-1);
            }   
            if (strcmp("eu", argv[i]) == 0 || strcmp("EU", argv[i]) == 0) {
                distanceType = EU;
            } else if (strcmp("js", argv[i]) == 0 || strcmp("JS", argv[i]) == 0) {
                distanceType = JS;            
            } else if (strcmp("ev", argv[i]) == 0 || strcmp("EV", argv[i]) == 0) {
              	 distanceType = EV;
	}
        }
        // unknown flag
        else {            
            printHelp();
            cerr << "Unknown flag: " << argv[i] << endl;
            exit(-1);
        }
    }
}

unsigned char* readData(const char * const filename, unsigned long long& n) {
	
	FILE *file;
	if (!(file = fopen(filename, "r"))) {
		printf("Unable to open file %s \n", filename);
		return NULL;
	}
	fseek (file , 0 , SEEK_END);
  	n = ftell (file);
	unsigned char *result = new unsigned char[n*2];
	rewind(file);
	if (n > fread(result, sizeof(unsigned char), n, file)) {
		printf("Error reading file %s \n", filename);
		fclose(file);
		delete[] result;
		return NULL;
	}
	fclose(file);
	return result;
}

int main(int argc, char *argv[]){
	parseParameters(argc,  argv);
	char *inputFile = argv[argc - 1];
	unsigned long long n = 0;
	unsigned char *str = readData(inputFile, n);
	if (str == NULL)
		return 0;
	vector<string> patternSet;
	if(patternFile.length()==0){
		variance* variance_obj;
        int* laenge = new int[2];
        laenge[0]=laenge[1]=dontcare+weight;
        variance_obj = new variance(number,laenge, weight);
        variance_obj->Silent(); 
        variance_obj->Improve(500);                            //5000 mal versuchen ein Pattern zu verbessern
        patternSet = variance_obj->GetPattern();
	}
	else
	parsePattern(patternFile, patternSet);
    double seq=0;
    int length=0;
    double dna=0;
    char c;
    for(unsigned long long i=0; i<n;i++){
        c=toupper(str[i]);
        if(c=='>'){
            while(str[++i]!='\n');
                seq++;
        }
        if(isalpha(c)){
            switch(c){
                case 'A': 
                    dna++;
                    break; 
                case 'C': 
                    dna++;
                    break;    
                case 'G': 
                    dna++;
                    break;              
                case 'T': 
                    dna++;
                    break; 
                case 'N': 
                    dna++;
                    break;
            }
            length++;
        }
    }
    if(dna/length>0.9){
        if(seq>0x7FFF){
            cout << "Too many sequences" << endl;
            exit(0);
        }
        if(weight>32){
            cout << "Maximum weight for DNA sequences is 32" << endl;
            exit(0);
        }
        cout << seq << " DNA sequences read of total length "<< length << endl;
    	if(n>=0xffffffff)
    		spacedDNA<unsigned long long>(patternSet, str,n, distanceType, threads, weight, dontcare, revComp, output);
    	else
    		spacedDNA<unsigned int>(patternSet, str,n, distanceType, threads,weight, dontcare, revComp, output);
    }
    else{
        if(seq>0xFFFF){
            cout << "Too many sequences" << endl;
            exit(0);
        }
         if(weight>12){
            cout << "Maximum weight for protein sequences is 12" << endl;
            exit(0);
        }
        cout << seq << " protein sequences read of total length "<< length << endl;
        if(distanceType==EV){
            cout << "Evolutionary distance not available for protein sequences"<< endl;
            exit(0);
        }
        if(n>=0xffffffff)
            spacedProt<unsigned long long>(patternSet, str,n, distanceType, threads, weight, dontcare, output);
        else
            spacedProt<unsigned int>(patternSet, str,n, distanceType, threads, weight, dontcare, output);
    }
}


