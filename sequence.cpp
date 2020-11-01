#include <string>
#include <map>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include "sequence.h"

using namespace std;

// Extracts the sequence for a text file.
string Sequence::getSequence(string fileName) {
    fstream file;
    string sequence;
    string line;
    file.open(fileName);
    while (getline(file, line)) {
        sequence.append(line);
    }
    return sequence;
}

// Constructor
Sequence::Sequence(string sequence, int kMerLength) {
    this -> inputSequence = sequence;
    this -> kMerLength = kMerLength;
    printf("kMer length = %i", kMerLength);
}

// Prints and returns the length of the input sequence
int Sequence::getLengthOfSequence() {
    printf("Sequence length = %i", int(inputSequence.size()));
    return inputSequence.size();
}

// Cuts the sequence at a specified length.
string Sequence::cutOffEndOfSequence() {
    printf("Cutting sequence at length = %i", maxSequenceLength);
    return inputSequence.substr(0, maxSequenceLength);
}

// This is used in the coming method to count the number of each nucleotide.
int Sequence::countNucleotide(int i) {
    vector<char> nucleotides = {'A', 'a', 'T', 't', 'C', 'c', 'G', 'g', 'N', 'n'};
    float upperCount = count(inputSequence.begin(), inputSequence.end(), nucleotides[2 * i]);
    float lowerCount = count(inputSequence.begin(), inputSequence.end(), nucleotides[(2 * i) + 1]);
    return upperCount + lowerCount;
}

// Gets the frequency of each nucleotide. withN says if N is to be counted.
vector<float> Sequence::getNucleotideFrequency(bool withN = false) {
    int repeats;
    if (withN == true) {
        repeats = 5;
    } else {
        repeats = 4;
    }
    for (int i = 0; i < repeats; i++) {
        int counter = countNucleotide(i);
        nucleotideFrequency.push_back(counter);
        nucleotideFrequency[i] /= inputSequence.size();
    }
    return nucleotideFrequency; // Returning a vector<float> with frequencies {'A', 'T', 'C', 'G'}
}

// This method will find the frequency of each dinucleotide in the sequence.
// It will return an unordered map of the dinucleotides with the frequency attached as the value.
// Clean the sequence before you attempt to get the dinucleotide frequencies unless you want N to be included.
unordered_map<string, float> Sequence::getDinucleotideFrequencyUMap() {
    string currentDinucleotide;
    for (int i = 0; i + 1 < int(inputSequence.size()); i++) {
        currentDinucleotide = inputSequence.substr(i, 2);
        if (dinucleotideFrequencyUMap[currentDinucleotide]) {
            dinucleotideFrequencyUMap[currentDinucleotide]++;
        } else {
            dinucleotideFrequencyUMap[currentDinucleotide] = 1;
        }
    }
    for (auto dinucleotide : dinucleotideFrequencyUMap) {
        float currentFrequency = dinucleotideFrequencyUMap[dinucleotide.first] / inputSequence.size();
        dinucleotideFrequencyUMap[dinucleotide.first] = currentFrequency;
    }
    return dinucleotideFrequencyUMap;
}

// This method breaks the sequence into small pieces and expresses the number of times that piece is found in the sequence.
unordered_map<string, vector<int>> Sequence::kMerizeWithCount() {
    const int inputSequenceSize = inputSequence.size();
    int i = 0;
    while (i + kMerLength <= inputSequenceSize) {
        string subSequence = inputSequence.substr(i, kMerLength);
        if (kMerMapCount.find(subSequence) == kMerMapCount.end()) {
            kMerMapCount.insert(pair<string, vector<int>>(subSequence, {1}));
        } else {
            kMerMapCount[subSequence][0]++;
        }
        i++;
    }
    return kMerMapCount;
}

// This method breaks the sequence into small pieces and tells you the starting position of each one of those pieces in the sequence.
unordered_map<string, vector<int>> Sequence::kMerizeWithPosition() {
    int i = 0;
    while (i + kMerLength <= int(inputSequence.size())) {
        string subSequence = inputSequence.substr(i, kMerLength);
        if (kMerMapPosition.find(subSequence) == kMerMapPosition.end()) {
            kMerMapPosition.insert(pair<string, vector<int>>(subSequence, {i}));
        } else {
            kMerMapPosition[subSequence].push_back(i);
        }
        i++;
    }
    return kMerMapPosition;
}

// Prints the kMerMap with count to the console.
void Sequence::printKmerMapCount() {
    for (auto keyValuePair : kMerMapCount) {
        cout << keyValuePair.first << " : " << keyValuePair.second[0] << endl;
    }
}

// Prints the kMerMap with position to the console.
void Sequence::printKmerMapPosition() {
    for (auto keyValuePair : kMerMapPosition) {
        for (int i = 0; i < int(keyValuePair.second.size()); i++) {
            cout << keyValuePair.first << " : " << keyValuePair.second[i] << endl;
        }
    }
}

// Prints the input sequence.
void Sequence::printSequence() {
    cout << inputSequence << endl;
}

// Turns each letter to upper case and screens out all non-nucleotide letters.
string Sequence::cleanSequence() {
    string sequence;
    for (int i = 0; i < int(inputSequence.size()); i++) {
        char currentChar = toupper(inputSequence[i]);
        if (currentChar == 'A' || currentChar == 'T' || currentChar == 'C' || currentChar == 'G') {
            sequence.push_back(inputSequence[i]);
        }
    }
    inputSequence = sequence;
    printf("Size of the clean sequence = %i", int(sequence.size()));
    return sequence;
}

// Creates a random nucleotide sequence where each nucleotide has a .25 chance of being picked.
string Sequence::createRandomSequenceDefault(int length) {
    int randomNumber;
    string sequence;
    int i = 0;
    while (i < length) {
        randomNumber = rand() % 4;
        if (randomNumber == 0) {
            sequence.push_back('A');
        } else if (randomNumber == 1) {
            sequence.push_back('T');
        } else if (randomNumber == 2) {
            sequence.push_back('C');
        } else if (randomNumber == 3) {
            sequence.push_back('G');
        }
    }
    printf("Size the random sequence = %i", int(sequence.size()));
    return sequence;
}

// Created a random nucleotide sequence using nucleotide frequencies.
string Sequence::createRandomSequenceWNF(int length) {
    int randomNumber;
    string sequence;
    vector<int> probabilityNumbers;
    for (int i = 0; i < 4; i++) {
        probabilityNumbers.push_back(nucleotideFrequency[i] * 100);
    }
    for (int i = 1; i < int(probabilityNumbers.size()); i++) {
        probabilityNumbers[i] += probabilityNumbers[i - 1];
    }
    for (int i = 0; i < length; i++) {
        randomNumber = rand() % 100;
        if (randomNumber < probabilityNumbers[0]) {
            sequence.push_back('A');
        } else if (randomNumber < probabilityNumbers[1] && randomNumber >= probabilityNumbers[0]) {
            sequence.push_back('T');
        } else if (randomNumber < probabilityNumbers[2] && randomNumber >= probabilityNumbers[1]) {
            sequence.push_back('C');
        } else {
            sequence.push_back('G');
        }
    }
    printf("Size the random sequence = %i", int(sequence.size()));
    return sequence;
}

// Creates a random nucleotide sequence using dinucleotide frequencies.
string Sequence::createRandomSequenceWDNF(int length) {
    string sequence;
    vector<float> probabilityNumbers;
    vector<string> dinucleotideCombinations;
    for (auto dinucleotide : dinucleotideFrequencyUMap) {
        probabilityNumbers.push_back(dinucleotide.second * 100);
        dinucleotideCombinations.push_back(dinucleotide.first);
    }
    for (int i = 1; i < int(probabilityNumbers.size()); i++) {
        probabilityNumbers[i] += probabilityNumbers[i - 1];
    }
    for (int i = 0; i  < length / 2; i++) {
        int randomNumber = rand() % 100;
        if (randomNumber < probabilityNumbers[0]) {
            sequence.append(dinucleotideCombinations[0]);
        } else {
            for (int j = 1; j < int(probabilityNumbers.size()); j++) {
                if (randomNumber > probabilityNumbers[j - 1] && randomNumber <= probabilityNumbers[j]) {
                        sequence.append(dinucleotideCombinations[j]);
                }
            }
        }
    }
    printf("Size the random sequence = %i", int(sequence.size()));
    return sequence;
}

// Given two kMerMaps, this method finds all the keys in one map that are also in the other map.
// withCount says if the kMerMap will be a count map or position map.
unordered_map<string, vector<int>> Sequence::getMatchList(unordered_map<string, vector<int>> kMerMap2, bool withCount) {
    unordered_map<string, vector<int>> matchUMap;
    unordered_map<string, vector<int>>kMerMap1;
    if (withCount == true) {
        kMerMap1 = kMerMapCount;
    } else {
        kMerMap1 = kMerMapPosition;
    }
    for(auto kMer2 : kMerMap2) {
        string currentSequence = kMer2.first;
        vector<int> currentNumericalData = kMer2.second; // Either count data or position data.
        if (kMerMap1.find(currentSequence) == kMerMap1.end()) {

        } else { // If the kMer is present in the other map
            if (matchUMap.find(currentSequence) == matchUMap.end()) {
                matchUMap.insert(pair<string, vector<int>>(currentSequence, currentNumericalData));
            }
        }
    }
    printf("Number of matches = %i", int(matchUMap.size()));
    return matchUMap;
}

