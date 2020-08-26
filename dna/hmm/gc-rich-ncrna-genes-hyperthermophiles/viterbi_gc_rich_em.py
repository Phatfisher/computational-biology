from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Viterbi import *
import sys

def print_gc_rich(viterbi_path):
  regions = []
  start = -1
  if viterbi_path[0] == 'gc_rich':
    start = 0

  for i in range(len(viterbi_path)-1):
    if viterbi_path[i] == 'normal' and viterbi_path[i+1] == 'gc_rich':
      start = i+1
    if viterbi_path[i] == 'gc_rich' and viterbi_path[i+1] == 'normal':
      end = i
      regions.append( (start, end) )
      start = -1

  if start != -1:
    # region hasn't yet ended
    regions.append( (start, range(len(viterbi_path)-1)) )
  print( 'GC-rich:', regions )

def em_train_and_print_gc_rich(genome_gc_vs_at, prior_probs, transition_probs, emission_probs, em_iters):
  def log_dict(d):
    return {a:numpy.log(b) for a,b in d.items()}

  log_prior_probs = log_dict(prior_probs)
  log_transition_probs = {a:log_dict(b) for a,b in transition_probs.items()}
  log_emission_probs = {a:log_dict(b) for a,b in emission_probs.items()}

  v = ViterbiDiscreteData(log_prior_probs, log_transition_probs, genome_gc_vs_at, log_emission_probs)
  viterbi_path, viterbi_log_prob = v.compute_viterbi_path_and_log_joint()
  print( 'Log probability', viterbi_log_prob )
  print()

  for i in range(em_iters):
    prior_probs = v.estimate_prior_probs_from_viterbi_path(viterbi_path)
    transition_probs = v.estimate_transition_probs_from_viterbi_path(viterbi_path)
    emission_probs = v.estimate_emission_probs_from_viterbi_path(viterbi_path)

    # ensure gc_rich refers to the state that more prefers S = {G,C}
    if emission_probs['normal']['S'] > emission_probs['gc_rich']['S']:
      emission_probs['normal'], emission_probs['gc_rich'] = emission_probs['gc_rich'], emission_probs['normal']

    log_prior_probs = log_dict(prior_probs)
    log_transition_probs = {a:log_dict(b) for a,b in transition_probs.items()}
    log_emission_probs = {a:log_dict(b) for a,b in emission_probs.items()}

    print(prior_probs)
    print(transition_probs)
    print(emission_probs)

    v = ViterbiDiscreteData(log_prior_probs, log_transition_probs, genome_gc_vs_at, log_emission_probs)
    viterbi_path, viterbi_log_prob = v.compute_viterbi_path_and_log_joint()
    print( 'Log probability', viterbi_log_prob )
    print()

  print_gc_rich(viterbi_path)

def usage():
  print('usage: <genome_fasta> [em]')

def main(args):
  if len(args) == 2:
    genome_fname, em_iters = args
    em_iters = int(em_iters)

    for record in SeqIO.parse(genome_fname, "fasta", alphabet=IUPAC.ambiguous_dna):
      print('Chromosome ', record.name)
      # S = {G,C}, W = {A,T}
      genome_gc_vs_at = str(record.seq).replace('G','S').replace('C','S').replace('A','W').replace('T','W')
      
      # random params:
      p = numpy.random.uniform(0.5, 1.)
      prior_probs = {'normal':p, 'gc_rich':1-p}
  
      p = numpy.random.uniform(0.5, 1.)
      q = numpy.random.uniform(0.5, 1.)
      transition_probs = { 'normal':{'normal':p, 'gc_rich':1-p}, 'gc_rich':{'normal':1-q, 'gc_rich':q}}
  
      p = numpy.random.uniform(0.5, 1.)
      q = numpy.random.uniform(0.5, 1.)
      emission_probs = { 'normal':{'S':1-p,'W':p}, 'gc_rich':{'S':q,'W':1-q} }
  
      print(prior_probs)
      print(transition_probs)
      print(emission_probs)
  
      em_train_and_print_gc_rich(genome_gc_vs_at, prior_probs, transition_probs, emission_probs, em_iters)
  else:
    usage()
    sys.exit(1)
    
if __name__=='__main__':
  main(sys.argv[1:])
