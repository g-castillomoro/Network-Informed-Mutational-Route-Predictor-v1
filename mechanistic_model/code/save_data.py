# AUTHOR: GASTON LUCA CASTILLO MORO
# DATE: DECEMBER 23, 2023

# This script saves the generated data with canonical names

import pickle

def save_sampling_data(cond_probability_biofilm, allowed_mutations, mutant_final_dgc_concentrations, y_baseline, pathway_name):
  with open(str(pathway_name) + r"_cond_probability_biofilm.pkl", 'wb') as f:
    pickle.dump(cond_probability_biofilm, f)
  with open(str(pathway_name) + r"_allowed_mutations.pkl", 'wb') as f:
    pickle.dump(allowed_mutations, f)
  with open(str(pathway_name) + r"_mutant_final_dgc_concentrations.pkl", 'wb') as f:
    pickle.dump(mutant_final_dgc_concentrations, f)
  with open(str(pathway_name) + r"_y_baseline.pkl", 'wb') as f:
    pickle.dump(y_baseline, f)

def save_final_probabilities(final_probabilities, pd_range, pe_range, pathway_name):
  with open(str(pathway_name) + r"_final_probabilities.pkl", 'wb') as f:
    pickle.dump(final_probabilities, f)
  with open(str(pathway_name) + r"_pd_range.pkl", 'wb') as f:
    pickle.dump(p_range, f)
  with open(str(pathway_name) + r"_pd_range.pkl", 'wb') as f:
    pickle.dump(p_range, f)
