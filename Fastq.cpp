#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include "Fastq.h"

using namespace std;

// Constructor to create a fastq object with the opened file
Fastq::Fastq(string fileName) {
    fastqFile.open(fileName);
    cout << "Using file: " << fileName << endl;
}

// Destructor
Fastq::~Fastq() {
    cout << "Destructor" << endl;
}

// Extracts the entire sequence from the fastq file.
string Fastq::getSequence() {
    string line;
    int idCounter = 1;
    int phredCounter = 0;
    int sequenceLineCounter = 0;
    int aLineCounter = 0;
    // Goes through the entire file, line by line and if the line contains a sequence, that line is pushed into the sequence.
    // The counters you we are to ensure that line is indeed a sequence and not some other line.
    while (getline(fastqFile, line)) {
        if (count(line.begin(), line.end(), '@') == 1 && aLineCounter == sequenceLineCounter) {
            idCounter = 1;
            phredCounter = 0;
            sequenceLineCounter = 0;
            aLineCounter = 0;
        } else if (line[0] == 'A' || line[0] == 'a' || line[0] == 'T' || line[0] == 't' || line[0] == 'C' || line[0] == 'c' || line[0] == 'G' || line[0] == 'g' && idCounter == 1 && phredCounter == 0) {
            sequence.append(line);
            sequenceLineCounter++;
        } else if (line[0] == '+') {
            phredCounter++;
        } else {
            aLineCounter++;
        }
    }
    return sequence;
}

// Gets the sequence with the start position of each kmer in each block of sequence.
unordered_map<string, vector<int>> Fastq::getSequenceWithKmerPosition(int kMerLength) {
    unordered_map<string, vector<int>> positionUM;
    string currentSequence;
    vector<int> currentPositionVector;
    string line;
    // These counters are needed for keeping track of position within the fastq while reading it.
    int idCounter = 0;
    int phredCounter = 0;
    int sequenceLineCounter = 0;
    int aLineCounter = 0;
    // This goes through the file, line by line and extracts the sequence like above.
    // When done with that sequence block, this goes through that sequence block and breaks it into kmers.
    // Starting position of each kmer is recorded.
    while (getline(fastqFile, line)) {
        if (count(line.begin(), line.end(), '@') == 1 && aLineCounter == sequenceLineCounter) {
            currentSequence.clear();
            currentPositionVector.clear();
            idCounter = 1;
            phredCounter = 0;
            sequenceLineCounter = 0;
            aLineCounter = 0;
        } else if (line[0] == 'A' || line[0] == 'a' || line[0] == 'T' || line[0] == 't' || line[0] == 'C' || line[0] == 'c' || line[0] == 'G' || line[0] == 'g' && idCounter == 1 && phredCounter == 0) {
            currentSequence.append(line);
            sequenceLineCounter++;
        } else if (line[0] == '+') {
            phredCounter++;
        } else if (phredCounter == 1) {
            // Breaks the current sequence into kmers and plugs the positions that the kmer is found in the position unordered map.
            int i = 0;
            while (i + kMerLength <= int(currentSequence.size())) {
                string subSequence = currentSequence.substr(i, kMerLength);
                if (positionUM.find(subSequence) == positionUM.end()) {
                    positionUM.insert(pair<string, vector<int>>(subSequence, {i}));
                    i++;
                } else {
                    positionUM[subSequence].push_back(i);
                    i++;
                }
            }
            aLineCounter++;
        }
    }
    // This is a unordered map of each kmer and a vector<int> saying the position of each kmer start.
    // For example, if kmer A is found to start at position 4, 29, and 306 in the sequences then
    // in the unordered map we would be kmer A and a vector {4,29,306}. Note that the starting positions don't have to be
    // in the same sequence block.
    return positionUM;
}

// Closes the fastq file.
void Fastq::closeFastq() {
    fastqFile.close();
}

