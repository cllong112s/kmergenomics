#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

class Fastq {
public:
    fstream fastqFile;
    string fileName;
    string sequence;

    // Open a fastq with the name fileName
    Fastq(string fileName);

    // Destructor
    ~Fastq();

    // Extract the sequence from the fastq
    string getSequence();

    // Goes through each sequence block and records the position of each kmer in the sequence
    unordered_map<string, vector<int>> getSequenceWithKmerPosition(int kMerLength);

    // Creates a fastq file with a certain number of blocks, each sequence being a certain length with a certain dinucleotide frequency
    void createRandomFastq(string fileName, int numberOfSequences, int sequenceLength, vector<int> sequenceDinucleotideFrequency);

    // Closes the fastq that we opened in the constructor
    void closeFastq();
};
