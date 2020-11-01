#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <fstream>

using namespace std;

class Sequence {
public:
    string inputSequence; // The sequence we want to create an object around.
    int sequenceLength; // The length of the input sequence.
    int kMerLength; // The length of kMer you wish to use.
    int maxSequenceLength; // The max length you wish the sequence to be. You have the option to cut the sequence down if you wish.
    unordered_map<string, vector<int>> matchUMap; // The list of kMers in both sequence A and sequence B.
    vector<float> nucleotideFrequency; // A vector of the nucleotide frequencies. {'A', 'T', 'C', G', 'N' (optional)}
    unordered_map<string, float> dinucleotideFrequencyUMap; // An unordered map of the dinucleotides and their frequencies in the given sequence.
    unordered_map<string, vector<int>> kMerMapCount; // An unordered map of the kMers found in the sequence and the number of times each kMer was found.
    unordered_map<string, vector<int>> kMerMapPosition; // An unordered map of the kMers found in the sequence and the starting position at which it was found.
    // Allows the user to input a raw sequence from a text file
    string getSequence(string);
    // Creates the sequence object
    Sequence(string, const int);
    // Gets the length of the input sequence.
    int getLengthOfSequence();
    // Trims the sequence to a certain length.
    string cutOffEndOfSequence();
    // Counts how many times each nucleotide occurs. This is used inside of a for loop.
    int countNucleotide(int i);
    // Counts each different nucleotides, then does count / sequenceLength.
    vector<float> getNucleotideFrequency(bool withN);
    // Counts each different dinucleotide combination, then does / (sequenceLength / 2).
    unordered_map<string, float> getDinucleotideFrequencyUMap();
    // Breaks down the sequence into smaller sub-sequences with count using a sliding window.
    unordered_map<string, vector<int>> kMerizeWithCount();
    // Breaks down the sequence into small sub-sequences with position of first kMer position of the kMer using a sliding window.
    unordered_map<string, vector<int>> kMerizeWithPosition();
    // Prints kMer : Count on each line.
    void printKmerMapCount();
    // Prints kMer : Starting position 1, starting position 2,...
    void printKmerMapPosition();
    // Prints the sequence.
    void printSequence();
    // Takes all the non-sequence characters out, ex) Takes all N's out.
    string cleanSequence();
    // Creates a random sequence using equal nucleotide frequencies.
    string createRandomSequenceDefault(int length);
    // Creates a random sequence using nucleotide frequencies.
    string createRandomSequenceWNF(int length);
    // Creates a random sequence using dinucleotide frequencies.
    string createRandomSequenceWDNF(int length);
    // Produces an unordered map of the kMers found in both sequences with their numerical data attached in a vector. withCount says if you are using a kMerMapCount or kMerMapPosition
    unordered_map<string, vector<int>> getMatchList(unordered_map<string, vector<int>> kMerMap2, bool withCount);
};
