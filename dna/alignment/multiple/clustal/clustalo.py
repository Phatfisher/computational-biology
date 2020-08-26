# Note: you will need to install clustal-omega (a package for your
# operating system, not for python) to run this.
# E.g., on Fedora linux, run
# sudo dnf install clustal-omega

from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import sys

def run_clustalo(fasta_fname):
  clustalo_cline = ClustalOmegaCommandline("clustalo", infile=fasta_fname, outfile='clustalo_output.aln', force=True, seqtype='DNA')
  # execute clustal:
  clustalo_cline()
  
  # read in the clustal output:
  alignment = AlignIO.read(open('clustalo_output.aln'), 'fasta')

  # print the alignment:
  print(alignment)
  print()

  print('Consensus sequence:')
  summary_alignment = AlignInfo.SummaryInfo(alignment)
  consensus = summary_alignment.gap_consensus()
  print(consensus)
  
def main(args):
  if len(args) == 1:
    fasta_fname = args[0]

    run_clustalo(fasta_fname)
  else:
    print('clustalo <fasta>')

if __name__ == '__main__':
  main(sys.argv[1:])
