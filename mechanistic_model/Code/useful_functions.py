# AUTHOR: GASTON LUCA CASTILLO MORO
# DATE: DECEMBER 23, 2023

import numpy as np
from scipy.integrate import solve_ivp
from numpy import random
import pickle

# First, we define the classes for each pathway

class castillom_wsp:
  def __init__(self):
    self.name = 'wsp'
    self.index_dgc = 4
    self.num_rxns = 6
    self.num_reactants = 9

    # here we generate a Gene Exclusion Matrix (GEEM), an num_rxns x num_rxns matrix, where the entry in row i and column j, 
    # g_ij, is 1 iff reaction i and reaction j cannot change simultaneously because they would involve mutations in two 
    # separate genetic components. Otherwise g_ij = 0. This will be used for sampling later

    self.geem_matrix = [[0, 0, 0, 0, 1, 1],
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 1, 1],
                        [0, 0, 0, 0, 0, 0],
                        [1, 1, 1, 0, 0, 0],
                        [1, 0, 1, 0, 0, 0]]

  # function that defines the ODE system for Wsp
  def odes(self, t, y, reaction_rates, S, C_m):

    # we define the concentrations of each species from y values
    ABD = y[0]
    ABD_m = y[1]
    ABD_mp = y[2]
    E_p = y[3]
    R_p = y[4]
    E = y[5]
    F = y[6]
    F_p = y[7]
    R = y[8]

    # we retrieve reaction rates to simplify syntax
    r1 = reaction_rates[0]
    r2 = reaction_rates[1]
    r3 = reaction_rates[2]
    r4 = reaction_rates[3]
    r5 = reaction_rates[4]
    r6 = reaction_rates[5]

    # we compute the derivatives, denoting the derivative as y_prime.
    # the indexes correspond to the previous variable assignment from y.
    y_prime = np.zeros(len(y))

    y_prime[0] = r2*F_p*ABD_m-r1*ABD*C_m # dABD/dt
    y_prime[1] = r1*ABD*C_m-r2*F_p*ABD_m-r3*S*ABD_m+r4*E*ABD_mp # dABD_m/dt
    y_prime[2] = -r4*E*ABD_mp + r3*S*ABD_m # dABD_mp/dt
    y_prime[3] = r4*E*ABD_mp - r5*E_p*R -r6*E_p*F # dE_p/dt
    y_prime[4] = r5*R*E_p # dR_p/dt
    y_prime[5] = r6*E_p*F-r4*E*ABD_mp+r5*R*E_p # dE/dt
    y_prime[6] = -r6*E_p*F +r2*F_p*ABD_m # dF/dt
    y_prime[7] = r6*E_p*F -r2*F_p*ABD_m # dF_p/dt
    y_prime[8] = -r5*R*E_p-.01*R # dR/dt

    return y_prime # we return the computed derivatives for ODE solver

class castillom_yfi:
  def __init__(self):
    self.name = 'yfi'
    self.index_dgc = 6
    self.num_rxns = 4
    self.num_reactants = 7

    # here we generate a Gene Exclusion Matrix (GEEM), an num_rxns x num_rxns matrix, where the entry in row i and column j, 
    # g_ij, is 1 iff reaction i and reaction j cannot change simultaneously because they would involve mutations in two 
    # separate genetic components. Otherwise g_ij = 0. This will be used for sampling later

    self.geem_matrix = [[0, 0, 1, 1],
                        [0, 0, 0, 1],
                        [1, 0, 0, 0],
                        [1, 1, 0, 0]]

    self.p_i = 0 # this is probability that an indel occurs in the pathway's genes

  # function that defines the ODE system for Yfi
  def odes(self, t, y, reaction_rates, S, C_m):

    # we define the concentrations of each species from y values
    R = y[0]
    NR = y[1]
    B = y[2]
    B_a = y[3]
    BR = y[4]
    N = y[5]
    NN = y[6]

    # we retrieve reaction rates to simplify syntax
    r1 = reaction_rates[0]
    r2 = reaction_rates[1]
    r3 = reaction_rates[2]
    r4 = reaction_rates[3]

    # we compute the derivatives, denoting the derivative as y_prime.
    # the indexes correspond to the previous variable assignment from y.
    y_prime = np.zeros(len(y))

    y_prime[0] = -r2*R*B_a-r3*R*N # dR/dt
    y_prime[1] = r3*R*N # dNR/dt
    y_prime[2] = -r1*S*B # dB/dt
    y_prime[3] = r1*S*B-r2*R*B_a # dB_a/dt
    y_prime[4] = r2*B_a*R # dBR/dt
    y_prime[5] = -r3*R*N-r4*N*N # dN/dt
    y_prime[6] = r4*N*N # dNN/dt

    return y_prime # we return the computed derivatives for ODE solver

class castillom_morA:
  def __init__(self):
    self.name = 'morA'
    self.index_dgc = 1
    self.num_rxns = 3
    self.num_reactants = 6

    # here we generate a Gene Exclusion Matrix (GEEM), an num_rxns x num_rxns matrix, where the entry in row i and column j, 
    # g_ij, is 1 iff reaction i and reaction j cannot change simultaneously because they would involve mutations in two 
    # separate genetic components. Otherwise g_ij = 0. This will be used for sampling later

    self.geem_matrix = [[0, 0, 0],
                        [0, 0, 0],
                        [0, 0, 0]]

    self.p_i = 0 # this is probability that an indel occurs in the pathway's genes

  # function that defines the ODE system for morA
  def odes(self, t, y, reaction_rates, S, C_m):

    # we define the concentrations of each species from y values
    D = y[0]
    D_a = y[1]
    P = y[2]
    PD = y[3]
    C = y[4]
    DC = y[5]

    # we retrieve reaction rates to simplify syntax
    r1 = reaction_rates[0]
    r2 = reaction_rates[1]
    r3 = reaction_rates[2]

    # we compute the derivatives, denoting the derivative as y_prime.
    # the indexes correspond to the previous variable assignment from y.
    y_prime = np.zeros(len(y))

    y_prime[0] = -r1*D*S; # dD/dt
    y_prime[1] = r1*D*S; # dD_a/dt
    y_prime[2] = -r2*P*D_a; # dP/dt
    y_prime[3] = r2*P*D_a# dPD/dt
    y_prime[4] = -r3*C*D_a; # dC/dt
    y_prime[5] = r3*C*D; # dDC/dt

    return y_prime # we return the computed derivatives for ODE solver

class castillom_sia:
  def __init__(self):
    self.name = 'sia'
    self.index_dgc = 6
    self.num_rxns = 4
    self.num_reactants = 7

    # here we generate a Gene Exclusion Matrix (GEEM), an num_rxns x num_rxns matrix, where the entry in row i and column j, 
    # g_ij, is 1 iff reaction i and reaction j cannot change simultaneously because they would involve mutations in two 
    # separate genetic components. Otherwise g_ij = 0. This will be used for sampling later

    self.geem_matrix = [[0, 0, 1, 1],
                        [0, 0, 0, 0],
                        [1, 0, 0, 0],
                        [1, 0, 0, 0]]

    self.p_i = 0 # this is probability that an indel occurs in the pathway's genes

  # function that defines the ODE system for Sia
  def odes(self, t, y, reaction_rates, S, C_m):

    # we define the concentrations of each species from y values
    A = y[0]
    A_a = y[1]
    C_p = y[2]
    C = y[3]
    B = y[4]
    D = y[5]
    DC = y[6]

    # we retrieve reaction rates to simplify syntax
    r1 = reaction_rates[0]
    r2 = reaction_rates[1]
    r3 = reaction_rates[2]
    r4 = reaction_rates[3]

    # we compute the derivatives, denoting the derivative as y_prime.
    # the indexes correspond to the previous variable assignment from y.
    y_prime = np.zeros(len(y))

    y_prime[0] = -r1*A*S # dA/dt
    y_prime[1] = r1*A*S-r2*A_a*C_p # dA_a/dt
    y_prime[2] = r3*C*B - r2*A_a*C_p # dC_p/dt
    y_prime[3] = -r3*C*B + r2*A_a*C_p - r4*D*C # dC/dt
    y_prime[4] = -r3*C*B # dB/dt
    y_prime[5] = -r4*D*C # dD/dt
    y_prime[6] = r4*D*C # dDC/dt

    return y_prime # we return the computed derivatives for ODE solver

# We now proceed with the Bayesian Sampling protocol and the required functions needed to calculate the conditional probability of biofilm phenotype

def bayesian_sampling_protocol(total_runs, pathway):
  run_number = 0
  y_baseline = np.zeros(total_runs)
  max_num_mutations = 3**(pathway.num_rxns) # we establish upper bound for number of possible mutations
  allowed_mutations = return_allowed_mutations(pathway, max_num_mutations)
  mutant_final_dgc_concentrations = np.zeros((len(allowed_mutations), total_runs)) # this will constitute a matrix later on

  for run_number in range(0, total_runs):
    # sample reaction rates from [0.1, 100]
    reaction_rates = [10**(4*random.rand()-2) for x in
                      range(0, pathway.num_rxns)]
    r_rates_save = reaction_rates.copy() # save reaction rates for later
    # sample initial conditions from [0, 10]. Range chosen based on previous literature
    initial_conditions = [10*random.rand() for x in
                          range (0, pathway.num_reactants)]
    # sample initial condition for signal S from [0, 10]
    S = 10*random.rand()
    # initialize final time for use later in solving with random factors
    final_time = 0
    if pathway.name == 'wsp': # if wsp pathway, we also need to sample C_m
      C_m = 10*random.rand()
    else:
      C_m = 0 # if not wsp pathway, we assign 0 and it will be ignored later

    # baseline ODE solver

    tf = 1.0 # initial final time for ODE solver
    distance = 1.0 # initialized distance between last two solution points to determine convergence to steady-state of solution
    tolerance = 10**(-8)

    while distance > tolerance:
      result_ode = solve_ivp(pathway.odes, [0, tf], initial_conditions, method = 'LSODA', args=(reaction_rates, S, C_m)) # result_ode is a dictionary
      distance = sum([(result_ode.y[i][-1] - result_ode.y[i][-2])**2 for i in range(0, pathway.num_reactants)]) # squared Euclidean distance between last two solution points
      y_baseline[run_number] = result_ode.y[pathway.index_dgc][-1]
      final_time = result_ode.t[-1]
      tf *= 2

    y_mutated = mutation_sampling(pathway, r_rates_save, initial_conditions, final_time, S, C_m) # allowed_mutations and y_mutated is built so that the ith entry of y_mutated
                                                                                                 # corresponds to the one produced by the mutation in the ith entry in allowed_mutations
    for i in range(0, len(y_mutated)):
      mutant_final_dgc_concentrations[i][run_number] = y_mutated[i] # this way we have the same numeration for allowed mutations and their corresponding DGC concentration

  cond_probability_biofilm = np.zeros(len(allowed_mutations)) # we initialize array to save conditional probabilities for each mutation
  for i in range(0, len(allowed_mutations)):
    biofilm_phenotype = 0 # counter for the times that mutated pathway yields higher concentration of DGC
    for j in range(0, total_runs):
      if mutant_final_dgc_concentrations[i][j] > y_baseline[j]:
        biofilm_phenotype += 1
    cond_probability_biofilm[i] = biofilm_phenotype / total_runs # this yields the conditional probability that a biofilm phenotype arises given a particular mutation in the pathway

  return cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline

def return_allowed_mutations(pathway, max_num_mutations):
  allowed_mutations = []
  for i in range (0, max_num_mutations):
    # decomposing into base-3 representation is a computationally-efficient way to generate all possible combinations of n = num_rxns of {-1, 0, 1} (even though here we use {0, 1, 2} for convenience)
    mutation_combinations = np.array(base3_representation(i, pathway.num_rxns))
    if check_combination(pathway, mutation_combinations):
      allowed_mutations.append(i) # this constitutes an array with the decimal representation of the mutation combination
  return allowed_mutations
  
def mutation_sampling(pathway, r_rates_save, initial_conditions, tf, S, C_m):
  max_num_mutations = 3**(pathway.num_rxns) # we establish upper bound for number of possible mutations
  y_steady_state_mutated = []
  allowed_mutations = return_allowed_mutations(pathway, max_num_mutations)

  for i in range (0, len(allowed_mutations)):
      mutation = base3_representation(allowed_mutations[i], pathway.num_rxns)
      factors = [10**((-2)*random.rand()), 1, 10**(2*random.rand())] # randomly sample factors for disabling mutations between [0.01, 1], nothing, and for enabling mutations between [1, 100] (in that order)
      reaction_rates = [(r_rates_save[i] * factors[int(mutation[i])]) for i in range(0, pathway.num_rxns)] # here we alter baseline reaction rates accordingly
      result_ode = solve_ivp(pathway.odes, [0, tf], initial_conditions, method = 'LSODA', args=(reaction_rates, S, C_m)) # result_ode is a dictionary
      y_steady_state_mutated.append(result_ode.y[pathway.index_dgc][-1])

  return y_steady_state_mutated

def check_combination(pathway, mutation_combinations):
  # we check whether the combination is feasible with our criteria using the GEEM matrix of the pathway
  for i in range (0, pathway.num_rxns):
    for j in range (0, pathway.num_rxns): # we start from i + 1 as matrix is symmetric and g_ii is always 0
      if pathway.geem_matrix[i][j] == 1:
        if int(mutation_combinations[i]) != 1 and int(mutation_combinations[j]) != 1: # there is an 'illegal' combination as they have both mutated
          return False
        else: # there is at least one that is not 1, so it is feasible
          pass
  return True # the mutation is feasible under our criteria

def base3_representation(number, num_rxns):
  n_base3 = np.zeros(num_rxns)
  for j in range(1, num_rxns + 1):
      number, n_base3[j-1] = divmod(number, 3)
  return n_base3[::-1] # this last step returns a reversed array

# Next, is the function to calculate the final probabilities of each pathway

def calculate_final_mutation_probabilities(pathway, p_d, p_e, cond_prob_biofilm, allowed_mutations):
  pathway_probability = 0
  for i in range(0, len(allowed_mutations)):
    mutation = base3_representation(allowed_mutations[i], pathway.num_rxns) - np.array([1]*pathway.num_rxns)
    mutation_probability = 1 # we initialize the probability that a mutation occurs in pathway
    for j in range(0, pathway.num_rxns): # this for loop calculates the probability of mutation occuring
      if mutation[j] > 0:
        mutation_probability *= p_e
      elif mutation[j] < 0:
        mutation_probability *= p_d
      else:
        mutation_probability *= (1 - p_d - p_e)
    pathway_probability += (cond_prob_biofilm[i] * mutation_probability) # this adds the contributing probability of each mutation to that of the whole pathway
  return pathway_probability
