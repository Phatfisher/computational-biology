from needleman_wunsch import *

def safe_slice(list_or_none, start=None, end=None, gap=None):
  if list_or_none is None:
    return None
  else:
    return list_or_none[start:end:gap]

def advance_column_ending_here_rightward(A, B, prev_column, col_index, upper_seed_row=None):
  col_size = len(A)+1

  next_column = numpy.zeros(col_size) - numpy.inf
  if upper_seed_row is not None:
    next_column[0] = upper_seed_row[col_index+1]

  for i in range(col_size):
    next_column[i] = max(next_column[i], prev_column[i] + gap_penalty)
    if i > 0:
      next_column[i] = max(next_column[i], next_column[i-1] + gap_penalty)
      next_column[i] = max(next_column[i], prev_column[i-1] + score_pair(A[i-1], B[col_index]) )

  return next_column

def advance_column_starting_here_leftward(A, B, next_column, col_index, lower_seed_row=None):
  return advance_column_ending_here_rightward(A[::-1], B[::-1], next_column[::-1], len(B)-col_index-1, upper_seed_row=safe_slice(lower_seed_row,gap=-1))[::-1]

def get_column_ending_here(A, B, n, upper_seed_row=None, left_seed_col=None):
  assert(left_seed_col is None or len(A)+1 == len(left_seed_col) )
  assert(upper_seed_row is None or len(B)+1 == len(upper_seed_row) )

  assert(n >=0 and n<len(B)+1)

  col_size = len(A)+1

  column = left_seed_col
  if left_seed_col is None:
    column = numpy.arange(0, col_size)*gap_penalty

  if n == 0:
    return column

  for col in range(n):
    column = advance_column_ending_here_rightward(A, B, column, col, upper_seed_row)

  return column

def get_column_starting_here(A, B, n, lower_seed_row=None, right_seed_col=None):
  assert(right_seed_col is None or len(A)+1 == len(right_seed_col) )
  assert(lower_seed_row is None or len(B)+1 == len(lower_seed_row) )

  return get_column_ending_here(A[::-1], B[::-1], len(B)-n, upper_seed_row=safe_slice(lower_seed_row, gap=-1), left_seed_col=safe_slice(right_seed_col, gap=-1))[::-1]

def get_column_passing_through_here(A, B, n, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None):
  return get_column_ending_here(A, B, n, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col) + get_column_starting_here(A, B, n, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

def get_row_ending_here(A, B, n, upper_seed_row=None, left_seed_col=None):
  assert(left_seed_col is None or len(A)+1 == len(left_seed_col) )
  assert(upper_seed_row is None or len(B)+1 == len(upper_seed_row) )

  return get_column_ending_here(B, A, n, upper_seed_row=left_seed_col, left_seed_col=upper_seed_row)

def get_row_starting_here(A, B, n, lower_seed_row=None, right_seed_col=None):
  assert(right_seed_col is None or len(A)+1 == len(right_seed_col) )
  assert(lower_seed_row is None or len(B)+1 == len(lower_seed_row) )

  return get_column_starting_here(B, A, n, lower_seed_row=right_seed_col, right_seed_col=lower_seed_row)

def get_row_passing_through_here(A, B, n, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None):
  return get_row_ending_here(A, B, n, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col) + get_row_starting_here(A, B, n, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

def needleman_wunsch_linear_space_helper(A, B, result, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None, top_row_index=0, left_column_index=0):
  middle_column_index = len(B)//2
  middle_column = get_column_passing_through_here(A, B, middle_column_index, upper_seed_row, left_seed_col, lower_seed_row, right_seed_col)

  # row cell at which middle column has best score passing through:
  best_val, best_row_index = max([(b,a) for a,b in enumerate(middle_column) ])

  corresponding_row = get_row_passing_through_here(A, B, best_row_index, upper_seed_row, left_seed_col, lower_seed_row, right_seed_col)

  if len(A) <= 2 or len(B) <= 2:
    base_case_indices = [ (i+top_row_index,j+left_column_index) for i,j in nw_backtrack(A, B, len(A), len(B)) ]
    result.extend(base_case_indices)
  else:
    # fixme: put into result

    # recurse on upper left:
    top_row = get_row_starting_here(A, B, 0, lower_seed_row, right_seed_col)[:middle_column_index+1]
    left_col = get_column_starting_here(A, B, 0, lower_seed_row, right_seed_col)[:best_row_index+1]
    bottom_row = get_row_ending_here(A, B, best_row_index, upper_seed_row, left_seed_col)[:middle_column_index+1]
    right_col = get_column_ending_here(A, B, middle_column_index, upper_seed_row, left_seed_col)[:best_row_index+1]
    needleman_wunsch_linear_space_helper(A[:best_row_index], B[:middle_column_index], result, top_row, left_col, bottom_row, right_col, top_row_index, left_column_index)

    # recurse on lower right:
    top_row = get_row_starting_here(A, B, best_row_index, lower_seed_row, right_seed_col)[middle_column_index+1:]
    left_col = get_column_starting_here(A, B, middle_column_index, lower_seed_row, right_seed_col)[best_row_index+1:]
    bottom_row = get_row_ending_here(A, B, len(A), upper_seed_row, left_seed_col)[middle_column_index+1:]
    right_col = get_column_ending_here(A, B, len(B), upper_seed_row, left_seed_col)[best_row_index+1:]
    needleman_wunsch_linear_space_helper(A[best_row_index+1:], B[middle_column_index+1:], result, top_row, left_col, bottom_row, right_col, top_row_index+best_row_index+1, left_column_index+middle_column_index+1)

def needleman_wunsch_linear_space(A, B):
  result = []
  needleman_wunsch_linear_space_helper(A, B, result)
  return sorted(result)

def main(args):
  if len(args) == 2:
    seq_a, seq_b = args

    print( needleman_wunsch_linear_space(seq_a, seq_b) )
  else:
    print('needleman_wunsch_linear_space <seq_a> <seq_b>')

if __name__ == '__main__':
  main(sys.argv[1:])
