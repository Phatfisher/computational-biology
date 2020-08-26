#include "Clock.hpp"
#include <fstream>
#include <string>

using namespace serang_lab;

constexpr unsigned int NUCLEOTIDES_PER_BLOCK = sizeof(unsigned long) * 8 / 2;

char nucleotide_complement(char nuc) {
  if (nuc == 'G')
    return 'C';
  if (nuc == 'C')
    return 'G';
  if (nuc == 'A')
    return 'T';
  if (nuc == 'T')
    return 'A';
  return -1;
}

std::string complement(const std::string & genome) {
  std::string result;

  for (char c : genome)
    result += nucleotide_complement(c);

  return result;
}

int main() {
  // A file of 4E6 G A T and C characters (the contents are
  // unimportant):
  std::ifstream fin("bigger.txt");
  std::string genome;

  Clock c;

  char base;
  while (fin >> base)
    genome += base;

  std::string comp = complement(genome);

  c.ptock();

  for (unsigned int i=0; i<2*NUCLEOTIDES_PER_BLOCK; ++i)
    std::cout << genome[i] << " ";
  std::cout << "..." << std::endl;
  for (unsigned int i=0; i<2*NUCLEOTIDES_PER_BLOCK; ++i)
    std::cout << comp[i] << " ";
  std::cout << "..." << std::endl;

  return 0;
}
