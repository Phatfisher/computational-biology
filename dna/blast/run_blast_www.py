from Bio.Blast import NCBIWWW
import sys

def main(args):
  if len(args) == 1:
    dna_seq = args[0]

    result = NCBIWWW.qblast("blastn", "nt", dna_seq)
    print( result.read() )

  else:
    print('usage: blastn <dna_seq>')
    
if __name__=='__main__':
  main(sys.argv[1:])

