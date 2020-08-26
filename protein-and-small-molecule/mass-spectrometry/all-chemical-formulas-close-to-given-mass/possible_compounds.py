import numpy, sys

def make_name(element_name, copies):
  if copies == 0:
    return ''
  if copies == 1:
    return element_name
  else:
    return element_name + str(copies)

def print_possible_compounds(goal_mass, element_to_mass, epsilon):
  compound_to_mass = {'':0.0}
  # generate all compounds with masses <= goal_mass+epsilon:
  for e,m in element_to_mass.items():
    max_copies_of_element = int(numpy.ceil( (goal_mass+epsilon)/m ))
    possible_element_contributions = {make_name(e,copies):(copies*m) for copies in range(max_copies_of_element)}
    
    new_compound_to_mass = {}
    for name_a,mass_a in compound_to_mass.items():
      for name_b,mass_b in possible_element_contributions.items():
        new_mass = mass_a+mass_b
        if new_mass <= goal_mass + epsilon:
          new_compound_to_mass[name_a+name_b] = new_mass

    compound_to_mass = new_compound_to_mass

  # print compounds within epsilon of goal_mass:
  for compound,mass in compound_to_mass.items():
    if numpy.fabs(mass - goal_mass) < epsilon:
      print('{:<30} {:.6}'.format(compound, mass))

def main(args):
  if len(args) == 3:
    goal_mass, elements_fname, epsilon = args
    goal_mass = float(goal_mass)
    # fixme: float
    element_to_mass = {}
    for line in open(elements_fname).readlines():
      symbol, mass, element = line.split()
      element_to_mass[symbol] = float(mass)
    epsilon = float(epsilon)

    print_possible_compounds(goal_mass, element_to_mass, epsilon)

  else:
    print('usage: possible_compounds <goal_mass_of_compound> <elements_filename> <mass_epsilon>')
    
if __name__=='__main__':
  main(sys.argv[1:])
