/**
 * Sebastian Porras
 *
 * An implementation of a progressive sequence
 * alignment for protein sequences.
 */

#include "msa.h"
#include "matrix.h"

using namespace std;
    
int main(int argc, char **argv){

    if (argc == 1) {
        cout << "Provide a fasta file" << endl;
        return CLI_ERROR;
    }

    vector<Sequence> seqs = read_fasta_file(argv[FILENAME]);

    auto StartTimeRef = chrono::high_resolution_clock::now();

    //convert blosum into a direct access array 
    vector<int> subMatrix = make_sub_matrix(); 

    auto make_sub_matrix = chrono::high_resolution_clock::now();
    float time1 = chrono::duration_cast<chrono::nanoseconds>(make_sub_matrix - StartTimeRef).count();
    float time_fin1 = 1e-9 * time1;
    
    auto calc_dist_start = chrono::high_resolution_clock::now();

    //will hold all distances between sequence pairs 
    vector<float> distanceMatrix = vector<float>(seqs.size() * seqs.size());  
    
    MPI_Init(&argc, &argv);

    //calculate similarity matrix between all pairs of sequences 
    calc_distances(seqs.size(), seqs, subMatrix, distanceMatrix);

    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {

        auto calc_dist_end = chrono::high_resolution_clock::now();
        float calc_dist_end_time = chrono::duration_cast<chrono::nanoseconds>(calc_dist_end - calc_dist_start).count();
        float calc_dist_end_real = 1e-9 * calc_dist_end_time;
        
        // Assign each sequence to its own cluster
        vector<vector<Sequence>> clusters;
        for (int i = 0; i < (int) seqs.size(); ++i) {
            vector<Sequence> singleCluster(1, seqs[i]);
            clusters.push_back(singleCluster);
        }
    
        auto upgma_start = chrono::high_resolution_clock::now();
   
        UPGMA(clusters, distanceMatrix, subMatrix);

        auto FinishTimeRef = chrono::high_resolution_clock::now();
        float upgmaRef = chrono::duration_cast<chrono::nanoseconds>(FinishTimeRef - upgma_start).count();
        float upgmaTime = 1e-9 * upgmaRef;
    
        float TotalTimeRef = chrono::duration_cast<chrono::nanoseconds>(FinishTimeRef - StartTimeRef).count();
        float time = 1e-9 * TotalTimeRef;

        cout << "Load_BLOSUM() seconds: " << fixed << time_fin1 << 
        setprecision(9) << "\n"; 

        cout << "Create_Matrix (s): " << fixed << calc_dist_end_real << 
        setprecision(9) << "\n"; 
    
        cout << "UPGMA (s): " << fixed << upgmaTime << 
        setprecision(9) << "\n"; 
        
        cout << "Total (s): " << fixed << time << 
        setprecision(9) << "\n"; 
        
        cout << argv[FILENAME] << "\n"; 

        //print_seqs(clusters); 
    }

    MPI_Finalize();

    return 0; 
}

/**
 * Go through the blosum matrix and 
 * convert it into a direct access array 
 * based on ASCII characters of the different 
 * possible chemicals. 
 * 
 * Return (vector<int>):
 * The direct access array 
*/
vector<int> make_sub_matrix(void) {

    //the order of chemicals stored in the blosum matirx 
    //string aOrder= "ARNDCQEGHILKMFPSTWYV";
    int aOrder[20] = {65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87, 89, 86};
    vector<int> subMatrix(MATRIX_SIZE, 0);
    __m256i offset = _mm256_set1_epi32(ASCII_OFFSET);
    __m256i rowLen = _mm256_set1_epi32(ROW_LEN);


    for (int i = 0; i < NUM_LETTERS; i++) {
        for (int j = 0; j < NUM_LETTERS - 8; j += 8) {
            __m256i acidI = _mm256_set1_epi32(aOrder[i]); 
            __m256i matPosI = _mm256_add_epi32(acidI, offset);
            __m256i matPosRow = _mm256_mullo_epi32(matPosI, rowLen);

            __m256i acidJ = _mm256_loadu_si256((__m256i*) &aOrder[j]);
            __m256i matPosJ = _mm256_add_epi32(acidJ, offset);

            __m256i indices = _mm256_add_epi32(matPosRow, matPosJ);

            int indexArr[8];
            _mm256_storeu_si256((__m256i*)indexArr, indices);

            __m256i scores = _mm256_loadu_si256((__m256i*) &blosum[i][j]);
            int scoreArr[8];
            _mm256_storeu_si256((__m256i*)scoreArr, scores);

            subMatrix[indexArr[0]] = scoreArr[0];
            subMatrix[indexArr[1]] = scoreArr[1];
            subMatrix[indexArr[2]] = scoreArr[2];
            subMatrix[indexArr[3]] = scoreArr[3];
            subMatrix[indexArr[4]] = scoreArr[4];
            subMatrix[indexArr[5]] = scoreArr[5];
            subMatrix[indexArr[6]] = scoreArr[6];
            subMatrix[indexArr[7]] = scoreArr[7];
        }
       
        for (int j = 16; j < NUM_LETTERS; j++) {
            //take the ASCII value of the char and add the correct alignment 
            //score at this position based off the blosum matrix.
            subMatrix[(aOrder[i] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i][j]; 
        }    
    }
 
    return subMatrix; 
}

/**
 * Perform the UPGMA clustering. Each sequence starts out in its own cluster.
 * UPGMA which involves finding the two closest clusters, finding the
 * two most similar sequences and then aliging everything within the clusters.
 * 
 * Repeat this process until there is only one cluster left. 
 * 
 * clusters (vector<vector<Sequence>>&): 
 * distanceMatrix (vector<float>&): 
 * subMatrix (vector<int>&): 
*/
void UPGMA(vector<vector<Sequence>> &clusters, 
        vector<float>& distanceMatrix, vector<int>& subMatrix) {

    int numClusters = clusters.size(); // iterate until there is 1 cluster
    const int numSeqs = numClusters; //track how many points we can compare

    while (numClusters > 1) {

        vector<Sequence> cToMerge1; //pointer to the next cluster to merge
        int idxC1 = 0; //will store the position in the list of clusters
        vector<Sequence> cToMerge2; //pointer to the other cluster to merge
        int idxC2 = 0;//will store the 2nd position in the list of clusters
    
        //locate two closest clusters based on average linkage 
        find_closest_clusters(numClusters, clusters, numSeqs, distanceMatrix, 
                              cToMerge1, &idxC1, cToMerge2, &idxC2); 

        //find two closest sequences within the cluster and align 
        align_clusters(cToMerge1, cToMerge2, subMatrix);

        // check which idx is greater so order is not messed up when removing
        if (idxC1 > idxC2) {
            clusters.erase(clusters.begin() + idxC1);
            clusters.erase(clusters.begin() + idxC2);
        } else {
            clusters.erase(clusters.begin() + idxC2);
            clusters.erase(clusters.begin() + idxC1);
        }

        // collapse old clusters into one new cluster
        vector<Sequence> newCluster = merge_clusters(cToMerge1, cToMerge2);

        //add the merged cluster back to the pile 
        clusters.push_back(newCluster); 
        numClusters -= 1;
    }
}

/**
 * merge_clusters
 * _______________
 * 
 * Take two clusters known to be closest to one 
 * another and merge them together. 
 * cToMerge1 (vector<Sequence>&): First cluster to merge 
 * cToMerge2 (vector<Sequence>&): Second cluster to merge 
 * 
 * Return (vector<Sequence>):
 * The newly merged cluster
 * 
*/
vector<Sequence> merge_clusters(vector<Sequence>& cToMerge1,
                                vector<Sequence>& cToMerge2) {

    // collapse old clusters into new cluster
    vector<Sequence> newCluster;

    int i; //add each sequence to the new cluster 
    for (i = 0; i < (int) cToMerge1.size() - 3; i += 4) {
        newCluster.push_back(cToMerge1[i]);
        newCluster.push_back(cToMerge1[i + 1]);
        newCluster.push_back(cToMerge1[i + 2]);
        newCluster.push_back(cToMerge1[i + 3]);
    }

    //add any remaining sequences 
    for (; i < (int) cToMerge1.size(); ++i) {
        newCluster.push_back(cToMerge1[i]);
    }

    //repeat but for the second cluster 
    for (i = 0; i < (int) cToMerge2.size() - 3; i += 4) {
        newCluster.push_back(cToMerge2[i]);
        newCluster.push_back(cToMerge2[i + 1]);
        newCluster.push_back(cToMerge2[i + 2]);
        newCluster.push_back(cToMerge2[i + 3]);
    }

    for (; i < (int) cToMerge2.size(); ++i) {
            newCluster.push_back(cToMerge2[i]);
    }

    return newCluster; 
}


/**
 * find_closest_clusters
 * _____________________
 * 
 * Checks through all pairs of remaining clusters 
 * and finds the two that are closest to one another. 
 * 
 * Updates values which keep track of the indices of the 
 * two clusters that are most similar to each other. 
 * 
 * numClusters (int): clusters left to merge 
 * clusters (vector<vector<Sequence>>&): contains all clusters  
 * numSeqs (int): number of sequences to compare 
 * distanceMatrix (vector<float>&): holds all distances between sequences  
 * cToMerge1 (vector<Sequence>*): pointer to first cluster to merge 
 * idxC1 (int*): index of first cluster within the clusters vector 
 * cToMerge2 (vector<Sequence>*): pointer to second cluster to merge 
 * idxC2 (int*): index of second cluster within the clusters vector 
 * 
*/
void find_closest_clusters(int numClusters, vector<vector<Sequence>> &clusters,
                           int numSeqs, vector<float>& distanceMatrix, 
                           vector<Sequence>& cToMerge1, int* idxC1, 
                           vector<Sequence>& cToMerge2, int* idxC2) {


    //keep track of the current smallest distance 
    float mostSimilar = DBL_MAX; 
 
    //iterate through all pairs of clusters 
    #pragma omp parallel
    {

        int localC1Idx = -1;    
        int localC2Idx = -1;

        float localMostSimilar = DBL_MAX; 
 
        #pragma omp for 
        for (int i = 0; i < numClusters; ++i) {
            for (int j = i + 1; j < numClusters; ++j) {
                
                //measure how close the these pair of clusters are 
                float dist = mean_difference(clusters[i], clusters[j], 
                        numSeqs, distanceMatrix);
                
                //perform the reduction over this region 
                if (dist < localMostSimilar) {
                    localMostSimilar = dist; 
                    localC1Idx = i;    
                    localC2Idx = j;
                }
            }
        }

        //update the record of what two clusters are closest between threads
        #pragma omp critical
        {
        if (localMostSimilar < mostSimilar) {
            //record the new smallest distance to be compared 
            mostSimilar = localMostSimilar;

            //keep track of which clusters will need to be merged 
            cToMerge1 = clusters[localC1Idx];
            cToMerge2 = clusters[localC2Idx];

            //also keep track of their indices so they can be removed 
            *idxC1 = localC1Idx;
            *idxC2 = localC2Idx;
        }
        }
    }
}


/*
 * mean_difference
 * ________________
 * 
 Find the mean difference bewteen two clusters using UPGMA.
 *
 * The difference between two clusters is defined as the
 * average difference between each pair of points to every
 * other point. https://en.wikipedia.org/wiki/UPGMA
 * 
 * c1 (vector<Sequence>&): the first cluster
 * c2 (vector<Sequence>&): the second cluster 
 * numSeqs (int): The number of sequences in the distance array 
 * distanceMatrix (vector<float>): records the distances between all sequences
 * 
 * Return (float): 
 */
 
float mean_difference(vector<Sequence> &c1, vector<Sequence> &c2,
        const int numSeqs, vector<float>& distanceMatrix){

    float mean = 0.0; //will store average distance 
    const int c1Size = c1.size(); // record the size of each cluster
    const int c2Size = c2.size();
    int chunkCount = numSeqs / 8;

    // take each sequence in the cluster and add up the differences
    #pragma omp parallel for reduction(+:mean)
    for (int i = 0; i < c1Size; ++i) {

        Sequence seq1 = c1[i]; 
        int seq1Index = seq1.index; //will be used to index into dist matrix
        
        for (int j = 0; j < c2Size; ++j) {

            Sequence seq2 = c2[j]; //the second sequence to compare against 
            int seq2Index = seq2.index; // the index to look up in the matrix 

            mean += sqrt(seq_to_seq_distance(seq1Index, seq2Index, distanceMatrix, 
                             chunkCount,  numSeqs));    
        }
    }

    return mean / (c1Size * c2Size);
}

float seq_to_seq_distance(int seq1Index, int seq2Index, vector<float>& distanceMatrix, 
                            int chunkCount, int numSeqs) {

    __m256 dist = _mm256_set1_ps(0); //will hold the sum of all distances 

    //iterate through the distance matrix and compare similarity 
 
    for (int k = 0; k < chunkCount; k++) {

        int vec1Index = seq1Index * numSeqs + k * 8; 
        int vec2Index = seq2Index * numSeqs + k * 8; 

        __m256 vec1Dists = _mm256_loadu_ps(&distanceMatrix[vec1Index]);
        __m256 vec2Dists = _mm256_loadu_ps(&distanceMatrix[vec2Index]);  

        __m256 diff = _mm256_sub_ps(vec1Dists, vec2Dists); // get distance
        __m256 square = _mm256_mul_ps(diff, diff); //square diff

        dist = _mm256_add_ps(dist, square);
    }

    //add up all the elements in the vector to get the total 
    dist = _mm256_hadd_ps(dist, dist);
    dist = _mm256_hadd_ps(dist, dist);
    __m256 distV2 = _mm256_permute2f128_ps(dist, dist, 1);
    dist = _mm256_add_ps(dist, distV2);

    //convert vector back into array 
    float* distSum = (float*) &dist;

    //add up the remaining values in case our seq isn't divisible by 8
    int remaining = numSeqs % 8; 
    int startAt = 8 * chunkCount;

    for (int i = 0; i < remaining; i++) {

        float dist = distanceMatrix[seq1Index * numSeqs + (startAt + i)] 
                    - distanceMatrix[seq2Index * numSeqs + (startAt + i)];
        dist *= dist; 
        distSum[0] += dist;
    }

    return distSum[0];
}


/** 
 * read_fasta_file
 * ________________
 * 
 * takes a fasta file and load contents into Sequence structs.
 * Will return a vector containing all the sequences. Will
 * exit with a FILE_ERROR if the file is not valid.
 * 
 * fileName(string): name of fasta file
 */
vector<Sequence> read_fasta_file(string fileName) {
    
    ifstream file(fileName); // open file

    if (!file.is_open()) {
        cout << "File could not be opened" << endl;
        exit(FILE_ERROR);
    }

    vector<Sequence> seqs; //store sequences here 
    string line; //current line we're on 
    string currentSeq;
    string currentId;
    int seqCount = 0;
    
    Sequence newSeq;
    while (getline(file, line)) {

        //indicates new sequence 
        if (line[0] == '>') {

            // our current seq is done, save it and clear space for next one 
            if (!currentSeq.empty()) {        
                newSeq.seq = currentSeq; 
                newSeq.id = currentId;
                newSeq.index = seqCount;
                seqs.push_back(newSeq);
                seqCount++;

                currentSeq.clear(); // start next seq fresh
            }
            
            //our new ID is the current line 
            currentId = line;

        } else { // all other lines are sequence information 
            currentSeq += line;
        }
    }
    // save the last sequence
    newSeq.seq = currentSeq;
    newSeq.id = currentId;
    newSeq.index = seqCount;

    seqs.push_back(newSeq);
    file.close();

    return seqs;
}

/**
 * print_seqs
 * __________
 * Prints out sequence IDs followed by actual sequence 
 * as per the FASTA format
*/
void print_seqs(vector<vector<Sequence>> clusters) {
    for (int i = 0; i < (int) clusters.size(); ++i)
    {
        for (int j = 0; j < (int) clusters[i].size(); j++)
        {
            cout << clusters[i][j].id << "\n" << 
                clusters[i][j].seq << endl;
        }
    }
}
