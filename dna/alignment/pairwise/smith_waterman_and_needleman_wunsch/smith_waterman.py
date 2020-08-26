import numpy
import sys
sys.setrecursionlimit(20000)

sys.path.append('../../../../utility/Memoized.py')
from Memoized import *

gap_penalty = -4
mismatch_penalty = -2
match_bonus = 1

nucleotides = ['G', 'A', 'T', 'C']

match_scores = {}
for nuc in nucleotides:
  match_scores[nuc] = match_bonus

def score_pair(nuc_a, nuc_b):
  if nuc_a == '-' or nuc_b == '-':
    return gap_penalty
  if nuc_a != nuc_b:
    return mismatch_penalty
  return match_bonus

@Memoized
def get_best_score_ending_here(A, B, index_a, index_b):
  if index_a == 0 or index_b == 0:
    return 0

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]

  # always include case to start at this cell
  best_scores_ending_before_here = [0]

  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a-1, index_b) + score_pair(nuc_a, '-'))
  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a, index_b-1) + score_pair('-', nuc_b))
  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b))

  return max(best_scores_ending_before_here)

@Memoized
def get_best_score_starting_here(A, B, index_a, index_b):
  return get_best_score_ending_here(A[::-1], B[::-1], len(A)-index_a, len(B)-index_b)

@Memoized
def get_best_score_pasing_through(A, B, index_a, index_b):
  return get_best_score_starting_here(A, B, index_a, index_b) + get_best_score_ending_here(A, B, index_a, index_b)

def print_matrix(A, B, score):
  mat = numpy.zeros( (len(A)+1, len(B)+1) )
  for i in range(len(A)+1):
    for j in range(len(B)+1):
      mat[i,j] = score(A,B,i,j)
  print(mat)

# backtrack:
def backtrack(A, B, index_a, index_b, score_ending_here_function):
  score_ending_here = score_ending_here_function(A, B, index_a, index_b)

  my_path = [(index_a, index_b)]

  if index_a == 0 or index_b == 0:
    return my_path

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]
   
  # Note: beware of floating point arithmetic ==:
  if score_ending_here_function(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b-1, score_ending_here_function)
  elif score_ending_here_function(A, B, index_a-1, index_b) + score_pair(nuc_a, '-') == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b, score_ending_here_function)
  elif score_ending_here_function(A, B, index_a, index_b-1) + score_pair('-', nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a, index_b-1, score_ending_here_function)
  else: # best path began here (upper, upper-left, and left are not good):
    return my_path

def print_alignment_from_path(A, B, path):
  alignment_a = ''
  alignment_b = ''
  for i in range(len(path)-1):
    index_a, index_b = path[i]
    next_index_a, next_index_b = path[i+1]

    if (index_a + 1, index_b + 1) == (next_index_a, next_index_b):
      # diagonal move:
      alignment_a += A[next_index_a-1]
      alignment_b += B[next_index_b-1]
    elif index_a == next_index_a:
      # vertical move
      alignment_a += '-'
      alignment_b += B[next_index_b-1]
    else:
      # horizontal move
      alignment_a += A[next_index_a-1]
      alignment_b += '-'

  print('alignment')
  print(alignment_a)
  print(alignment_b)

def main(args):
  if len(args) == 2:
    seq_a, seq_b = args

    print('Starting here:')
    print_matrix(seq_a, seq_b, get_best_score_starting_here)
    print()
  
    print('Ending here:')
    print_matrix(seq_a, seq_b, get_best_score_ending_here)
    print()
  
    print('Passing through:')
    print_matrix(seq_a, seq_b, get_best_score_pasing_through)
    print()
  
    best_score, best_index = max( [ (get_best_score_pasing_through(seq_a, seq_b, i,j),(i,j)) for i in range(len(seq_a)+1) for j in range(len(seq_b)+1) ] )
  
    print('best score', best_score, 'from best index', best_index)
  
    best_path = backtrack(seq_a, seq_b, *best_index, get_best_score_ending_here)[::-1]
    print('best path', best_path)
  
    print_alignment_from_path(seq_a, seq_b, best_path)
  
  else:
    print('smith_waterman <seq_a> <seq_b>')

if __name__ == '__main__':
  main(sys.argv[1:])
