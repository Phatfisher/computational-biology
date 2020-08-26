import sys

# fixme: molecular evolution model through time
class DNASubstitutionModel:
  # define in derived class
  def probability_of_y_given_x(self, x, y):
    pass

def DNASubstitutionModelFelsenstein81(DNASubstitutionModel):
  def __init__(self, pi_g, pi_a, pi_t):
    pi_c = 1 - pi_g - pi_a - pi_t
    self._nucleotide_to_pi = {'G':pi_g, 'A':pi_a, 'T':pi_t, 'C':pi_c}

    self._beta = 1. / (1 - pi_g**2 - pi_a**2 - pi_t**2 - pi_c**2)

  # x evolves to y
  def probability_of_y_given_x(self, x, y, time):
    val = self._nucleotide_to_pi[y] * (1 - numpy.exp(-self._beta*time))
    if x != y:
      return val
    else:
      return numpy.exp(-self._beta*time) + val

def single_base_pruning_log_probability(tree_node, name_to_sequence, substitution_model, base_index):
  # fixme:
  pass

def multi_base_pruning_log_probability(tree_node, name_to_sequence, substitution_model):
  # all seqs have same length (due to gaps in global alignment). use an arbitrary one to decide on length
  n = len(list(name_to_sequence.values())[0])

  return sum([ single_base_pruning_log_probability(tree_node, name_to_sequence, substitution_model, i) for i in range(n) ])

def main(args):
  if len(args) == 2:
    fasta_mult_alignment_fname, newick_tree_fname = args

    substitution_model = DNASubstitutionModelFelsenstein81(0.3, 0.2, 0.2, 0.3)

    # read in the multiple alignment fasta (e.g., from clustal output):
    name_to_sequence = {record.name:record.seq for record in AlignIO.read(open(fasta_mult_alignment_fname), 'fasta') }

    # fixme: read newick tree (using dict keys as leaves)
    

    print( 'Log Pr =', multi_base_pruning_log_probability(root_node, name_to_sequence, substitution_model) )
  else:
    print('usage: digest <multiple_alignment_fasta_file> {ecori,bamhi,hindiii}')
    
if __name__=='__main__':
  main(sys.argv[1:])
