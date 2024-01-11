# AUTHOR: GASTON LUCA CASTILLO MORO
# DATE: DECEMBER 23, 2023

import useful_functions as uf
import save_data as sd

def generate_all_data(pathway, num_samples = 10000, step_size_pd_range = 0.05, step_size_pe_range = 0.05):
  # num_samples refers to the number of samples for the Bayesian Sampling protocol & step_size_pd_range,step_size_pe_range refers to the step size 
  # in the generation of the set used for p_d and p_e, respectively
  
  cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline = uf.bayesian_sampling_protocol(num_samples, pathway)
  pd_range = 10 ** np.arange(-7, -1, step_size_p_range)
  pe_range = 10 ** np.arange(-7, -1, step_size_p_range)
  final_probabilities = np.zeros((len(pd_range), len(pe_range)))
  pd_counter = 0
  pe_counter = 0
  for p_d in pd_range:
    pe_counter = 0
    for p_e in pe_range:
      final_probabilities[pd_counter][pe_counter] = uf.calculate_final_mutation_probabilities(pathway, p_d, p_e, cond_probability_biofilm, allowed_mutations)
      pe_counter += 1
    pd_counter += 1
    print('pd_counter' + str(pd_counter))
  sd.save_sampling_data(cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline, pathway.name)
  sd.save_final_probabilities(final_probabilities, pd_range, len(pd_range), pathway.name)

if __name__ == '__main__':
  generate_all_data()
