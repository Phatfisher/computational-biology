from Bio import SeqIO, motifs
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

import pylab as P

import sys

def get_upstream_prefixes_and_downstream_suffixes(gb_file, chromosome_to_sequence, prefix_suffix_size):
  prefix_result = []
  suffix_result = []
  for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank"):
    chromosome = gb_record.name
    chromosome_sequence = chromosome_to_sequence[chromosome]
    for feature in gb_record.features:
      locus = feature.location
      if feature.type == 'CDS':
        # note: these may not handle regions that bridge the gap from
        # the end to the start of circular chromosomes

        if locus.strand == 1:
          # sense strand:
          start = locus.start-1
          five_prime_prefix = chromosome_sequence[start-prefix_suffix_size:start]

          end = locus.end+1
          three_prime_suffix = chromosome_sequence[end:end+prefix_suffix_size]
        else:
          # sense strand:
          end = locus.end+1
          five_prime_prefix = chromosome_sequence[end:end+prefix_suffix_size].reverse_complement()

          start = locus.end-1
          three_prime_suffix = chromosome_sequence[start-prefix_suffix_size:start].reverse_complement()

        # only allow things that have length prefix_size (for simplicity)
        if len(five_prime_prefix) == prefix_suffix_size:
          prefix_result.append(five_prime_prefix.upper())
        if len(three_prime_suffix) == prefix_suffix_size:
          suffix_result.append(three_prime_suffix.upper())
  return prefix_result, suffix_result

def print_consensus_and_draw_sequence_heatmap(sequences):
  # allow for ambiguity codes
  m = motifs.create(sequences, alphabet=IUPAC.ambiguous_dna)
  print(m.consensus)
  d = m.counts.normalize()
  P.imshow([d['G'], d['A'], d['T'], d['C']], interpolation=None)
  P.gca().set_yticks(range(4))
  P.gca().set_yticklabels(['G','A','T','C'])
  P.colorbar(orientation='horizontal')

def main(args):
  if len(args) != 2:
    print('usage: <fasta_file> <genbank_file>')
  else:
    fa_fname, gb_fname = args

    chromosome_to_sequence = {}
    # allow for ambiguity codes to match print_consensus_and_draw_sequence_heatmap (alphabets should be consistent)
    for record in SeqIO.parse(fa_fname, "fasta", alphabet=IUPAC.ambiguous_dna):
      # remove version (after . symbol):
      record_id_without_version = record.id.split('.')[0]
    
      chromosome_to_sequence[record_id_without_version] = record.seq
    
    prefix_and_suffix_length = 150
    prefixes, suffixes = get_upstream_prefixes_and_downstream_suffixes(gb_fname, chromosome_to_sequence, prefix_and_suffix_length)

    P.figure(1)
    print('prefix consensus: ', end='')
    print_consensus_and_draw_sequence_heatmap(prefixes)
    P.xlabel('Base (CDS is to the right)')

    P.figure(2)
    print('suffix consensus: ', end='')
    print_consensus_and_draw_sequence_heatmap(suffixes)
    P.xlabel('Base (CDS is to the left))')

    P.show()

    # save the prefixes and suffixes in case later processing is wanted:
    prefixes = [ SeqRecord(prefixes[i],id='PREF_{}'.format(i)) for i in range(len(prefixes)) ]
    SeqIO.write(prefixes, open('prefixes.fasta','w'), format='fasta')

    suffixes = [ SeqRecord(suffixes[i],id='SUFF_{}'.format(i)) for i in range(len(suffixes)) ]
    SeqIO.write(suffixes, open('suffixes.fasta','w'), format='fasta')

if __name__=='__main__':
  main(sys.argv[1:])
