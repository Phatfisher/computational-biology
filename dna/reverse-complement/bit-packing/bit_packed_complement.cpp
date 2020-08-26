#include "Clock.hpp"
#include <fstream>
#include <string>
#include <bitset>

using namespace serang_lab;

// Alternate ordering makes it so that we can complement with the ~ operator:
// G -> 0 = 00
// A -> 1 = 01
// T -> 2 = 10
// C -> 3 = 11
// -1uc (i.e., 255) otherwise
const unsigned char nuc_to_code[] = {255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 1, 255, 3, 255, 255, 255, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};

constexpr unsigned int NUCLEOTIDES_PER_BLOCK = sizeof(unsigned long) * 8 / 2;

// Rewritten to modify reference variables so that both result and
// num_blocks are provided (a return value could only return 1 unless
// we wrap it in std::pair<>):
void bit_pack(unsigned long * & result, unsigned long & num_blocks, const std::string & genome) {
  num_blocks = genome.size() / NUCLEOTIDES_PER_BLOCK;
  // If it doesn't divide easily (* should be faster than %), add one
  // more block:
  if ( num_blocks * NUCLEOTIDES_PER_BLOCK != genome.size() )
    ++num_blocks;

  result = (unsigned long*)calloc(num_blocks, sizeof(unsigned long));

  unsigned long i=0;
  // Go through all but the final block:
  for (unsigned long block=0; block<num_blocks-1; ++block) {
    unsigned long next_block = 0;
    
    for (unsigned char j=0; j<NUCLEOTIDES_PER_BLOCK; ++j, ++i) {
      next_block <<= 2;
      
      next_block |= nuc_to_code[ genome[i] ];
    }
    result[block] = next_block;
  }

  // Perform final block separately, and only check boundary of genome
  // in final block:
  unsigned long next_block = 0;
  for (unsigned char j=0; j<NUCLEOTIDES_PER_BLOCK; ++j, ++i) {
    next_block <<= 2;
    
    // This is the final block: check to make sure the boundaries are met:
    if (i < genome.size())
      next_block |= nuc_to_code[ genome[i] ];
  }
  result[num_blocks-1] = next_block;
}

unsigned long* packed_complement(const unsigned long*packed_genome, unsigned long num_blocks) {
  unsigned long * result = (unsigned long*)calloc(num_blocks, sizeof(unsigned long));
  for (unsigned long block=0; block<num_blocks; ++block)
    // With our new ordering of the nucleotides, bit-wise complement
    // will perform nucleotide complement!
    result[block] = ~packed_genome[block];

  // The above will also flip any unused bits in the last block, but
  // those are to be ignored anyway.
  return result;
}

int main() {
  // A file of 4E6 G A T and C characters (the contents are
  // unimportant):
  std::ifstream fin("bigger.txt");
  std::string genome;

  char base;
  while (fin >> base)
    genome += base;

  unsigned long*packed;
  unsigned long num_blocks;
  bit_pack(packed, num_blocks, genome);

  Clock c;

  unsigned long*comp = packed_complement(packed, num_blocks);
  
  c.ptock();
  
  for (unsigned int i=0; i<2*NUCLEOTIDES_PER_BLOCK; ++i)
    std::cout << genome[i] << " ";
  std::cout << "..." << std::endl;

  for (unsigned int i=0; i<2; ++i)
    std::cout << std::bitset<8*sizeof(unsigned long)>(packed[i]);
  std::cout << "..." << std::endl;

  for (unsigned int i=0; i<2; ++i)
    std::cout << std::bitset<8*sizeof(unsigned long)>(comp[i]);
  std::cout << "..." << std::endl;

  char ordered_nucleoties[] = {'G','A','T','C'}; // Order matches code above
  for (unsigned int i=0; i<2; ++i) {
    unsigned long mask = -1ul & ~((-1ul)>>2); // 11000....
    for (unsigned char j=0; j<NUCLEOTIDES_PER_BLOCK; ++j) {
      unsigned long shifted_bit_pair = comp[i] & mask;
      std::cout << ordered_nucleoties[ shifted_bit_pair >> 2*(NUCLEOTIDES_PER_BLOCK-j-1) ] << " ";
      mask >>= 2;
    }
  }
  std::cout << "..." << std::endl;

  return 0;
}
