import sys
sys.path.append('../alignment/pairwise/smith_waterman_and_needleman_wunsch/')
from smith_waterman import *

import numpy
from scipy.optimize import curve_fit
from collections import Counter
import pylab as P

def random_seq(n):
  nucs = ['G','A','T','C']
  result = ''
  for i in range(n):
    result += nucs[ numpy.random.randint(0,4) ]
  return result

def simulate_histogram(number_iterations):
  scores = []
  length=100
  for iter in range(number_iterations):
    seq_a = random_seq(length)
    seq_b = random_seq(length)

    best_score = max([ get_best_score_ending_here(seq_a, seq_b, i, j) for i in range(length+1) for j in range(length+1) ])
    scores.append(best_score)

  score_to_count = Counter(scores)
  sorted_scores = sorted(set(scores))
  return sorted_scores, score_to_count

number_simulations = 1000

# with gaps allowed:
P.figure(1)
sorted_scores_gapped, score_to_count_gapped = simulate_histogram(number_simulations)
P.bar(sorted_scores_gapped, [ score_to_count_gapped[score] for score in sorted_scores_gapped ], label='Gapped', alpha=0.7)

# disallow gaps:
P.figure(2)
gap_penalty = -numpy.inf
sorted_scores_ungapped, score_to_count_ungapped = simulate_histogram(number_simulations)
x_data = sorted_scores_ungapped
y_data = numpy.array([ score_to_count_ungapped[score] for score in sorted_scores_ungapped ]) / float(number_simulations)
P.bar(x_data, y_data, label='Ungapped', alpha=0.7)

# Fit to a gumbel extreme-value distribution:
def gumbel_pdf(x, mu, beta):
  z = (x - mu) / float(beta)
  return 1.0/beta * numpy.exp(-(z + numpy.exp(-z)))

theta_opt, theta_cov = curve_fit(gumbel_pdf, x_data, y_data, bounds=([-20.0,0.0], [20.0,20]))
mu_opt, beta_opt = theta_opt
print( ('mu', 'beta'), '=', theta_opt )

x_min = min(sorted_scores_ungapped)
x_max = max(sorted_scores_ungapped)
x_interpol = numpy.linspace(x_min, x_max, 500)
P.plot(x_interpol, gumbel_pdf(x_interpol, mu_opt, beta_opt), label='Best fit')

P.yscale('log')
P.legend()
P.show()
