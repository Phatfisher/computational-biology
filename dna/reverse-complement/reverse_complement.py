def main(args):
  if len(args) == 1:
    # fixme: get dna_seq from args[0]

    # fixme: compute reverse complement of dna_seq

    print(rev_comp_dna_seq)
  else:
    print('needleman_wunsch <dna_seq>')

if __name__ == '__main__':
  main(sys.argv[1:])
