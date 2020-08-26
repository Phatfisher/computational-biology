from Bio import SeqIO
from Bio.Alphabet import IUPAC
import sys

class StringSliceLinear:
  def __init__(self, s, i, a_or_b):
    self._tup = (s, i, a_or_b)

  def get_start_index(self):
    return self._tup[1]

  def __lt__(self, rhs):
    x, i, x_a_or_b = self._tup
    y, j, y_a_or_b = rhs._tup
    return x[i:] < y[j:]

  def from_different_string_as(self, rhs):
    x, i, x_a_or_b = self._tup
    y, j, y_a_or_b = rhs._tup

    # use XOR to see if exactly one of these bools is true:
    return x_a_or_b ^ y_a_or_b

  def longest_shared_prefix_size(self, rhs):
    x, i, x_a_or_b = self._tup
    y, j, y_a_or_b = rhs._tup
    
    x = x[i:]
    y = y[j:]
    for k in range(min(len(x), len(y))):
      if x[k] != y[k]:
        break
    # k is the index of the first mismatch
    return k

class StringSliceCircular(StringSliceLinear):
  def __lt__(self, rhs):
    x, i, x_a_or_b = self._tup
    y, j, y_a_or_b = rhs._tup

    x_prefix = x[i:]
    y_prefix = y[j:]

    x_suffix = x[:i]
    y_suffix = y[:j]

    # Note: could be performed without making a copy
    return x_prefix + x_suffix < y_prefix + y_suffix

  def longest_shared_prefix_size(self, rhs):
    x, i, x_a_or_b = self._tup
    y, j, y_a_or_b = rhs._tup

    x_prefix = x[i:]
    y_prefix = y[j:]

    x_suffix = x[:i]
    y_suffix = y[:j]

    x_circ = x_prefix + x_suffix
    y_circ = y_prefix + y_suffix

    for k in range(min(len(x_circ), len(y_circ))):
      if x_circ[k] != y_circ[k]:
        break
    # k is the index of the first mismatch
    return k

def lcs_linear(A, B):
  print('Building suffixes...')
  # Note: these suffixes could simply be tuples of the form (A[i:],
  # True); however, python will use a memory-hungry sort, which can
  # easily exhaust your RAM. By using a class instead, we force python
  # to use comparison sort.
  suffixes = [ StringSliceLinear(A,i,True) for i in range(len(A)) ] + [ StringSliceLinear(B,i,False) for i in range(len(B)) ]
  print('Sorting suffixes...')
  suffixes.sort()

  print('Finding matching prefixes...')
  best_match_len = -1
  best_match = None

  # find suffixes from A and B adjacent (i.e., they share the longest prefix)
  for i in range(len(suffixes)-1):
    string_slice_x = suffixes[i]
    string_slice_y = suffixes[i+1]
    if string_slice_x.from_different_string_as(string_slice_y):
      # one is from A the other is from B:
      match_len = string_slice_x.longest_shared_prefix_size(string_slice_y)

      # j is the index of the first mismatch
      if match_len > best_match_len:
        best_match_len = match_len
        start_index = string_slice_x.get_start_index()
        best_match = A[start_index:start_index+match_len]
        print('new best match', best_match)

  return best_match

def lcs_circular(A, B):
  print('Building suffixes...')
  suffixes = [ StringSliceCircular(A,i,True) for i in range(len(A)) ] + [ StringSliceCircular(B,i,False) for i in range(len(B)) ]
  print('Sorting suffixes...')
  suffixes.sort()

  print('Finding matching prefixes...')
  best_match_len = -1
  best_match = None

  # find suffixes from A and B adjacent (i.e., they share the longest prefix)
  for i in range(len(suffixes)-1):
    string_slice_x = suffixes[i]
    string_slice_y = suffixes[i+1]
    if string_slice_x.from_different_string_as(string_slice_y):
      # one is from A the other is from B:
      match_len = string_slice_x.longest_shared_prefix_size(string_slice_y)

      # j is the index of the first mismatch
      if match_len > best_match_len:
        best_match_len = match_len
        start_index = string_slice_x.get_start_index()
        best_match = A[start_index:start_index+match_len]
        if start_index + match_len > len(A):
          # best match bridges over end of A back to start of A:
          best_match += A[:start_index + match_len - len(A)]

        print('new best match', best_match)

  return best_match

def main(args):
  if len(args) == 3:
    genome_fname_a, genome_fname_b, lin_or_circ = args
    assert(lin_or_circ in ('linear', 'circular'))

    (genome_a,) = SeqIO.parse(genome_fname_a, "fasta", alphabet=IUPAC.ambiguous_dna)
    genome_a = str(genome_a.seq)
    (genome_b,) = SeqIO.parse(genome_fname_b, "fasta", alphabet=IUPAC.ambiguous_dna)
    genome_b = str(genome_b.seq)
    
    if lin_or_circ == 'linear':
      result = lcs_linear(genome_a, genome_b)
    else:
      result = lcs_circular(genome_a, genome_b)
    print(result)

  else:
    print('lcs <fasta_a> <fasta_b> {linear,circular}')
    
if __name__=='__main__':
  main(sys.argv[1:])
