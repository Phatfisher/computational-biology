import sys

def main(args):
  if len(args) == 1:
    # fixme: load genome from the FASTA file, whose filename is found
    # in args[0]

    # fixme: compute gc_percentage and at_percentage

    print( gc_percentage, at_percentage )

  else:
    print('usage: gc_content <fasta_file>')
    
if __name__=='__main__':
  main(sys.argv[1:])
