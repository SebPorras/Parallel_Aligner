#ifndef MSA_H
#define MSA_H

#include <string>
#include <vector> 
#include <map>
#include <memory>
#include <list> 
#include <cmath> 
#include <limits.h>
#include <iomanip>
#include <ios>
#include <iostream>
#include <float.h>
#include <fstream> 
#include <chrono>
using namespace std;

//export out blosum matrix 
extern int blosum[20][20];

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1; //for accessing argv 
const int NUM_LETTERS = 20; //20 amino acids 
const int ROW_LEN = 24; //ascii value of 'Y' - 'A'
const int MATRIX_SIZE = 600; // substition matrix 24 * 24 + 24
const int ASCII_OFFSET = -65; //offsets 'A' to have index of 0

//error codes 
const int CLI_ERROR = 1;
const int FILE_ERROR = 2;

const int GAP = -3; //penalty for adding a gap 

//A representation of a sequence 
struct Sequence {
    string seq; //actual sequence 
    string id; //sequence name 
    int index; //where the sequence is in the matrix  
};

float mean_difference(vector<Sequence>& c1, vector<Sequence>& c2, 
        const int numPoints, vector<float> distanceMatrix); 
vector<Sequence> read_fasta_file(string fileName); 
void UPGMA(vector<vector<Sequence>>& clusters, 
        vector<float>& distanceMatrix, vector<int>& subMatrix);
void print_seqs(vector<vector<Sequence>> clusters); 
vector<int> make_sub_matrix(void);
void find_closest_clusters(int numClusters, vector<vector<Sequence>> &clusters,
                           int numSeqs, vector<float>& distanceMatrix, 
                           vector<Sequence>& cToMerge1, int* idxC1, 
                           vector<Sequence>& cToMerge2, int* idxC2);
vector<Sequence> merge_clusters(vector<Sequence>& cToMerge1, 
                                vector<Sequence>& cToMerge2);
                                
#endif