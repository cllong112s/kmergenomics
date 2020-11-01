#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace std;

class Fasta {
public:
    fstream fastaFile;
    string fileName;
    string sequence;

    // Open a fasta with the name "fileName".
    Fasta(string fileName);

    // Destructor.
    ~Fasta();

    // Extracts the sequence from the fasta file.
    string getSequence();

    // Creates a fasta file with a random sequence.
    void createRandomFastaDefault(string, int, int);

    // Creates a Fasta file with a random sequence as defined by dinucleotide frequencies.
    void createRandomFastaWDNF(string, int, int, unordered_map<string, float>);

    // Closes the Fasta file.
    void closeFasta();

};

