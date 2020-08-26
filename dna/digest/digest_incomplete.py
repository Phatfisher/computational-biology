import sys

def main(args):
  if len(args) == 2:
    # fixme: load genome from the FASTA file, whose filename is found
    # in args[0]



    # fixme: using the enzyme provided in args[1], digest the genome
    # into a list of DNA fragments. you need not consider IUPAC
    # ambiguity codes.



    # print the list of DNA fragments.

    # Note: print in sorted order so that the output should exactly
    # match any correct implementation.
    for fragment in sorted(list_of_dna_fragments):
      print(fragment)

  else:
    print('usage: digest <fasta_file> {ecori,bamhi,hindiii}')
    
if __name__=='__main__':
  main(sys.argv[1:])
