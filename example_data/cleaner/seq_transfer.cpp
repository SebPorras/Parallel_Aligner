//seq_transfer.cpp

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
using namespace std; 

struct Sequence {
    std::string seq; //actual sequence 
    std::string id; //sequence name 
    int index; //where the sequence is in the matrix  
};


std::vector<Sequence> read_fasta_file(std::string fileName)
{

    std::vector<Sequence> seqs;

    std::ifstream file(fileName); // open file

    if (!file.is_open())
    {
        std::cout << "File could not be opened" << std::endl;
        exit(3);
    }

    std::string line;
    std::string currentSeq;
    std::string currentId;
    int seqCount = 0;
    Sequence newSeq;

    while (std::getline(file, line))
    {

        if (line[0] == '>')
        {
            if (!currentSeq.empty())
            { // save our seq
                newSeq.seq = currentSeq;
                newSeq.id = currentId;
                newSeq.index = seqCount;
                seqs.push_back(newSeq);
                seqCount++;

                currentSeq.clear(); // start next seq fresh
            }

            currentId = line;
        }
        else
        { // all other lines are sequences
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


int main(int argc, char **argv){

    if (argc == 1)
    {
        cout << "Provide a fasta file" << endl;
        return 2;
    }

    // data structure to hold our sequences
    vector<Sequence> seqs = read_fasta_file(argv[1]);
    string newFile = string(argv[2]) + string("_seqs_globin");
    ofstream output(newFile);

    for (int i = 0; i < atoi(argv[2]); ++i) {
        output << seqs[i].id << endl;
        output << seqs[i].seq << endl; 
    }

    output.close();

    return 0; 
}
