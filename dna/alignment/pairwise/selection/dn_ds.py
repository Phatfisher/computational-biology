from Bio.Seq import Seq
import sys

from Bio.codonalign.codonseq import CodonSeq,cal_dn_ds

# CCC ACU AUC GUU AAC GAU AGC UGG UCC UAC
# CCA ACA AUG GUU AAC GAC AGA UCG UCC UAU 
#   P   T I/M           D S/R W/S       Y

def hamming_distance(codon_a, codon_b):
  return sum([ 1 if codon_a[i] != codon_b[i] else 0 for i in range(3) ])

def is_same_amino_acid(codon_a, codon_b):
  #print(codon_a, Seq(codon_a).translate(), 'vs', codon_b, Seq(codon_b).translate())
  return Seq(codon_a).translate() == Seq(codon_b).translate()

def possible_synonymous_and_nonsynonymous_changes(codon_a):
  synon_changes = nonsynon_changes = 0
  for c in 'GATC':
    for i in range(3):
      # ignore case where we change a base to itself (i.e., no change)
      if codon_a[i] != c:
        codon_a_copy = list(codon_a)
        codon_a_copy[i] = c
        
        if is_same_amino_acid(codon_a, ''.join(codon_a_copy)):
          synon_changes += 1
        else:
          nonsynon_changes += 1
  return synon_changes, nonsynon_changes

def estimate_dN_dS_by_counting_synonymous_and_nonsynonymous_changes_sum_then_ratio(seq_a, seq_b):
  # seq_a and seq_b form an ungapped alignment

  total_synon_changes = total_nonsynon_changes = 0
  total_possible_synon_changes = total_possible_nonsynon_changes = 0
  
  # traverse in codons
  for i in range(0,len(seq_a),3):
    codon_a = seq_a[i:i+3]
    codon_b = seq_b[i:i+3]

    print(codon_a, codon_b)
    
    # check that codons are full length:
    if len(codon_a) == 3 and len(codon_b) == 3:

      diffs = hamming_distance(codon_a, codon_b)
      if diffs == 1:
        possible_synon, possible_nonsynon = possible_synonymous_and_nonsynonymous_changes(codon_a)
        print('possible', possible_synon, possible_nonsynon)
        total_possible_synon_changes += possible_synon
        total_possible_nonsynon_changes += possible_nonsynon

        if is_same_amino_acid(codon_a, codon_b):
          # synon
          print('synon')
          total_synon_changes += 1
        else:
          print('nonsynon')
          total_nonsynon_changes += 1

    print()

  return (total_nonsynon_changes / total_possible_nonsynon_changes) / (total_synon_changes / total_possible_synon_changes)

def estimate_dN_dS_by_counting_synonymous_and_nonsynonymous_changes_sum_weighted(seq_a, seq_b):
  # seq_a and seq_b form an ungapped alignment

  total_synon_changes = total_nonsynon_changes = 0
  total_possible_synon_changes = total_possible_nonsynon_changes = 0
  
  # traverse in codons
  for i in range(0,len(seq_a),3):
    codon_a = seq_a[i:i+3]
    codon_b = seq_b[i:i+3]

    print(codon_a, codon_b)
    
    # check that codons are full length:
    if len(codon_a) == 3 and len(codon_b) == 3:

      diffs = hamming_distance(codon_a, codon_b)
      if diffs == 1:
        possible_synon, possible_nonsynon = possible_synonymous_and_nonsynonymous_changes(codon_a)
        print('possible', possible_synon, possible_nonsynon)
        total_possible_synon_changes += possible_synon
        total_possible_nonsynon_changes += possible_nonsynon

        if is_same_amino_acid(codon_a, codon_b):
          # synon
          print('synon')
          total_synon_changes += 1 / (possible_synon / (possible_synon + possible_nonsynon))
        else:
          print('nonsynon')
          total_nonsynon_changes += 1 / (possible_nonsynon / (possible_synon + possible_nonsynon))

    print()

  return total_nonsynon_changes / total_synon_changes

def main(args):
  if len(args) == 2:
    seq_a, seq_b = args
    assert(len(seq_a) == len(seq_b))

    print( 'dN/dS (ratio of sums) ~', estimate_dN_dS_by_counting_synonymous_and_nonsynonymous_changes_sum_then_ratio(seq_a, seq_b) )
    print('------------------')
    print( 'dN/dS  (sum weighted) ~', estimate_dN_dS_by_counting_synonymous_and_nonsynonymous_changes_sum_weighted(seq_a, seq_b) )

    print('------------------')
    dN,dS=cal_dn_ds(CodonSeq(seq_a),CodonSeq(seq_b))
    print( 'dN/dS    (with model) ~', dN/dS )
    print()

  else:
    print('ka_ks <dna_start_seq> <dna_end_seq>')

if __name__ == '__main__':
  main(sys.argv[1:])
