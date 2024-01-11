import useful_functions as uf

def save_sampling_data(cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline, pathway_name):
  with open(str(pathway_name) + r"_cond_probability_biofilm.pkl", 'wb') as f:
    pickle.dump(cond_probability_biofilm, f)
  with open(str(pathway_name) + r"_allowed_mutations.pkl", 'wb') as f:
    pickle.dump(allowed_mutations, f)
  with open(str(pathway_name) + r"_mutant_final_dgc_concentrations.pkl", 'wb') as f:
    pickle.dump(mutant_final_dgc_concentrations, f)
  with open(str(pathway_name) + r"_y_baseline.pkl", 'wb') as f:
    pickle.dump(y_baseline, f)


def retrieve_sampling_data(pathway_name):
  with open(str(pathway_name) + r"_cond_probability_biofilm.pkl", 'rb') as f:
    cond_probability_biofilm = pickle.load(f)
  with open(str(pathway_name) + r"_allowed_mutations.pkl", 'rb') as f:
    allowed_mutations = pickle.load(f)
  with open(str(pathway_name) + r"_mutant_final_dgc_concentrations.pkl", 'rb') as f:
    mutant_final_dgc_concentrations = pickle.load(f)
  with open(str(pathway_name) + r"_y_baseline.pkl", 'rb') as f:
    y_baseline = pickle.load(f)

  return cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline


def generate_cond_prob_distributions(pathway, num_samples):
  cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline = bayesian_sampling_protocol(num_samples, pathway)
  save_sampling_data(cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline, pathway.name)


def save_final_probabilities(final_probabilities, p_range, length_p_range, pathway_name):
  with open(str(pathway_name) + r"_" + str(length_p_range) + r"length_final_probabilities.pkl", 'wb') as f:
    pickle.dump(final_probabilities, f)
  with open(str(pathway_name) + r"_" + str(length_p_range) + r"length_p_range.pkl", 'wb') as f:
    pickle.dump(p_range, f)


def retrieve_final_probabilities(pathway_name, length_p_range):
  path_final_prob = str(pathway_name) + r"_" + str(length_p_range) + r"length_final_probabilities.pkl"
  path_p_range = str(pathway_name) + r"_" + str(length_p_range) + r"length_p_range.pkl"
  slant_bar = r'/'
  drive_path = r"drive/MyDrive/LAB SMANIA/2023/Proyect__Probabilistic_model_of_pathway_usage_for_biofilm_phenotype_conversion/generated_data/pathway_mutation_probabilities_probability_space/"
  with open(drive_path + '/' + str(pathway_name) + '/' + path_final_prob , 'rb') as f:
    final_probabilities = pickle.load(f)
  with open(drive_path + '/' + str(pathway_name) + '/' + path_p_range, 'rb') as f:
    p_range = pickle.load(f)

  return final_probabilities, p_range


def generate_final_probabilities(pathway, cond_probability_biofilm, allowed_mutations, step_size_p_range):
  pd_range = 10 ** np.arange(-7, -1, step_size_p_range)
  pe_range = 10 ** np.arange(-7, -1, step_size_p_range)
  final_probabilities = np.zeros((len(pd_range), len(pe_range)))
  pd_counter = 0
  pe_counter = 0
  for p_d in pd_range:
    pe_counter = 0
    for p_e in pe_range:
      final_probabilities[pd_counter][pe_counter] = calculating_probabilities(pathway, p_d, p_e, cond_probability_biofilm, allowed_mutations)
      pe_counter += 1
    pd_counter += 1
    print('pd_counter' + str(pd_counter))
  save_final_probabilities(final_probabilities, pd_range, len(pd_range), pathway.name)


def calculating_probabilities(pathway, p_d, p_e, cond_prob_biofilm, allowed_mutations):
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

if __name__ == '__main__':
  heatmap_plotting_segmented_regions()
