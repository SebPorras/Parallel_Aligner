#ifndef MATRIX_H 
#define MATRIX_H

#include "msa.h"
#include <ostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <iterator>
#include <unordered_map>
using namespace std; 

float run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify, vector<int>& subMatrix);
float calculate_similarity(string seq1, string seq2); 
vector<int> create_matrix(string& seq1, string& seq2, 
        const int rows, const int cols, const size_t length, vector<int>& subMatrix);
void nw_seq_to_seq( string& seq1, string& seq2, string& aSeq1, 
        string& aSeq2, vector<int>& M, int rows, int cols);
vector<float> calc_distances(int numSeqs, vector<Sequence>& seqs, vector<int>& subMatrix);
void nw_on_group(vector<int>& M, int rows, int cols, 
        vector<Sequence>& group1, vector<Sequence>& group2);
void setup_group_alignment(vector<Sequence>& group1, 
        vector<Sequence>& group2, int g1Idx, int g2Idx, 
        vector<int>& subMatrix); 
void choose_seq_group_align(vector<Sequence>& group1, 
        vector<Sequence>& group2, vector<int>& subMatrix); 
void align_clusters(vector<Sequence>& cToMerge1, 
        vector<Sequence>& cToMerge2, vector<int>& subMatrix); 
#endif 