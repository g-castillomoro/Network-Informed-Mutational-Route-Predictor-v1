Software used
--------

- Python v3.10.12
- Numpy v1.23.5
- SciPy v1.11.4

Usage
--------

For users who want to **apply the model to any pathway of interest** have to create a **new Python class** in the provided code, with the same format as `castillom_wsp`, `castillom_sia`, etc. The necessary components are listed below. 

**Necessary attributes:**

- `self.name`
  - Name of the pathway
- `self.index_dgc`
  - Index in the y-vector (in the ODE method) that corresponds to DGC (or any other enzyme of interest)
- `self.num_rxns`
  - Number of reactions
- `self.num_reactants`
  - Number of reactants/species/nodes in the reaction network
- `self.geem_matrix`
  - As explained in code, this is a (num_rxns x num_rxns) matrix, such that g_ij = 1 iff reaction i and reaction j cannot change simultaneously because they would involve mutations in two separate genetic components. Otherwise g_ij = 0.


**Necessary method:**

- `odes(self, t, y, reaction_rates, extra_args)`
  - a method that defines the ODE & that is compatible with the ODE solver that is being used (in this case, `scipy.integrate.solve_ivp`). As it is dependent on the ODE solver used, it is advised to refer to the corresponding documentation for appropiate implementation (for `scipy.integrate.solve_ivp`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
  - `extra_args` encompasses quantities whose initial values (eg.: `S` or `C_m`) are important for the ODE system, but whose dynamics are not important for the problem and are thus not considered in the ODE system

Note
----
For the function arguments with fallback values assigned to them, the default values were used during the actual runs when generating the data
