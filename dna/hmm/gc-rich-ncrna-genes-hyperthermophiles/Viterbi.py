import numpy
from collections import Counter, defaultdict

class Viterbi:
  def __init__(self, log_prior_probs, log_transition_probs, data):
    self._log_prior_probs = log_prior_probs
    self._log_transition_probs = log_transition_probs

    self._states = list(log_prior_probs.keys())
    self._data = data

  def _log_prob_of_transition_from_state_s1_to_s2(self, s1, s2):
    return self._log_transition_probs[s1][s2]

  # define in derived classes
  def _log_prob_of_state_s_emitting_data_d(self, s, d):
    pass

  def compute_viterbi_path_and_log_joint(self):
    layers_from_left = []

    # first layer uses prior and emission:
    d = self._data[0]
    layer = {}
    for s in self._states:
      layer[s] = self._log_prior_probs[s] + self._log_prob_of_state_s_emitting_data_d(s,d)
  
    layers_from_left.append(layer)
    
    # layers 2, 3, ... use previous layer, emission, and all transitions:
    for i in range(len(self._data)-1):
      d = self._data[i+1]
  
      new_layer = {s:-numpy.inf for s in self._states}
  
      for start_state in self._states:
        for end_state in self._states:
          new_layer[end_state] = max( new_layer[end_state], layer[start_state] + self._log_prob_of_transition_from_state_s1_to_s2(start_state,end_state) + self._log_prob_of_state_s_emitting_data_d(end_state,d) )
  
      layer = new_layer
      layers_from_left.append(layer)
  
    viterbi_log_prob,viterbi_end_state = max([ (b,a) for a,b in layer.items()])
    
    # perform backtracking:
    viterbi_path = [viterbi_end_state]
    for i in range(len(self._data)-1)[::-1]:
      d = self._data[i+1]
      log_scores_and_start_statetates = []
      viterbi_end_state = max([ (layers_from_left[i][start_state] + self._log_prob_of_transition_from_state_s1_to_s2(start_state,viterbi_end_state),start_state) for start_state in self._states])[1]
      viterbi_path.append(viterbi_end_state)
  
    # path was built right-to-left; reverse so that it's left-to-right
    viterbi_path = viterbi_path[::-1]
    return viterbi_path, viterbi_log_prob

  def estimate_prior_probs_from_viterbi_path(self, viterbi_path):
    n = len(viterbi_path)
    counts = Counter(viterbi_path)
    result = {a:(b/n) for a,b in counts.items()}

    # put in all states in case some don't occur. start at prob=0
    for s in self._states:
      if s not in result:
        result[s] = 0.0

    return result

  def estimate_transition_probs_from_viterbi_path(self, viterbi_path):
    # there are n-1 transitions
    num_transitions = len(viterbi_path)-1

    transitions = [ (viterbi_path[i], viterbi_path[i+1]) for i in range(num_transitions) ]
    start_counts = Counter([a for a,b in transitions])
    transition_counts = Counter(transitions)

    result = defaultdict(dict)

    # put in all states in case some don't occur. start at prob=0
    for s1 in self._states:
      for s2 in self._states:
        result[s1][s2] = 0.0

    for start_state,end_state in transition_counts:
      # normalize so that they're conditional rather than joint probabilities
      result[start_state][end_state] = transition_counts[ (start_state,end_state) ] / start_counts[start_state]

    return result

  # define in derived class
  def estimate_log_emission_probs_from_viterbi_path(self, viterbi_path):
    pass

# when the data options are small enough for a table (e.g., nucleotides).
# for things like PhastCons scores, this may not be possible. In those
# cases, use a closed form likelihood for Pr(D_i|S_i=s_i).
class ViterbiDiscreteData(Viterbi):
  def __init__(self, log_prior_probs, log_transition_probs, data, log_emission_probs):
    Viterbi.__init__(self, log_prior_probs, log_transition_probs, data)
    self._log_emission_probs = log_emission_probs

  # override
  def _log_prob_of_state_s_emitting_data_d(self, s, d):
    return self._log_emission_probs[s][d]

  # override
  def estimate_emission_probs_from_viterbi_path(self, viterbi_path):
    n = len(viterbi_path)

    emissions = [ (viterbi_path[i], self._data[i]) for i in range(n) ]
    start_counts = Counter(viterbi_path)
    emission_counts = Counter(emissions)

    result = defaultdict(dict)

    # put in all states,data in case some don't occur. start at prob=0
    for s in self._states:
      for d in set(self._data):
        result[s][d] = 0.0

    for start_state,end_state in emission_counts:
      result[start_state][end_state] = emission_counts[ (start_state,end_state) ] / start_counts[start_state]

    return result
