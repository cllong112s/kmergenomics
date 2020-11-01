#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Fasta.h"

using namespace std;


// Opens a Fasta file and prints to the console what file name was used.
Fasta::Fasta(string fileName) {

    // Open the fasta file.
    fastaFile.open(fileName);

    // Prints to the console the file name used.
    cout << "Using file: " << fileName << endl;
}

// Destructor
Fasta::~Fasta() {
    cout << "Destructor" << endl;
}

// Extracts the sequence from the fasta file.
string Fasta::getSequence() {
    string line;
    // This goes through, line by line and asking 'Is there a '>' in this line?'.
    // If there isn't, then this is a sequence line and should be extracted as part of the sequence.
    while (getline(fastaFile, line)) {
        if (count(line.begin(), line.end(), '>') == 0) {
            sequence.append(line);
        }
    }
    // Returning the extracted sequence.
    return sequence;
}

// This creates a fasta file where each sequence is random, with .25 frequency for each nucleotide.
// >ID separates each sequence block.
void Fasta::createRandomFastaDefault(string fileName, int numberOfSequences, int sequenceLength) {
    ofstream randomFastaFile(fileName);
    // This goes through each sequence block.
    for (int i = 0; i < numberOfSequences; i++) {
        string randomSequence;
        // Here we are going through each piece of the sequence and picking a random number to assign a nucleotide.
        for (int j = 0; j < sequenceLength; j++) {
            int randomNumber = rand() % 3;
            if (randomNumber == 0) {
                randomSequence.push_back('A');
            } else if (randomNumber == 1) {
                randomSequence.push_back('T');
            } else if (randomNumber == 2) {
                randomSequence.push_back('C');
            } else if (randomNumber == 3) {
                randomSequence.push_back('G');
            }
        }
        // Here is where we put the >ID and sequence into the fasta file.
        randomFastaFile << ">ID" << endl;
        randomFastaFile << randomSequence << endl;
        randomSequence.clear();
    }
    randomFastaFile.close();
}

// This creates a fasta file with random blocks of sequence created with a dinucleotide frequency map.
// There are numberOfSequences blocks of sequenceLength length. >ID separates each block.
void Fasta::createRandomFastaWDNF(string fileName, int numberOfSequences, int sequenceLength, unordered_map<string, float> sequenceDinucleotideFrequencyUM) {
    // Here we are creating two vectors so we can pull dinucleotides by generating random numbers later.
    vector<float> dinucleotideFrequencies;
    vector<string> dinucleotides;
    float currentSum = 0;
    for (auto dinucleotide : sequenceDinucleotideFrequencyUM) {
        dinucleotides.push_back(dinucleotide.first);
        currentSum += dinucleotide.second * 100;
        cout << currentSum << endl;
        dinucleotideFrequencies.push_back(currentSum);
    }
    // Here we are going through each sequence block, creating a random sequence from the dinucleotide frequency unordered map and placing that
    // random sequence and >ID into the file.
    ofstream randomFastaFile(fileName);
    for (int i = 0; i < numberOfSequences; i++) {
        string randomSequence;
        // Here we need to implement a random sequence based on dinucleotide frequencies.
        for (int j = 0; j < sequenceLength / 2; j++) {
            int randomNumber = rand() % 100;
            for (int k = 1; k < int(dinucleotideFrequencies.size()); k++) {
                if (randomNumber <= dinucleotideFrequencies[k] && randomNumber > dinucleotideFrequencies[k - 1]) {
                    randomSequence.append(dinucleotides[k]);
                } else if (randomNumber <= dinucleotideFrequencies[k] && k == 1) {
                    randomSequence.append(dinucleotides[0]);
                }
            }

        }
        // Here is where we are putting the information into the fasta file.
        randomFastaFile << ">ID" << endl;
        randomFastaFile << randomSequence << endl;
        randomSequence.clear();
    }
    dinucleotides.clear();
    dinucleotideFrequencies.clear();
    randomFastaFile.close();
}

// This method closes the fasta file.
void Fasta::closeFasta() {
    fastaFile.close();
}


