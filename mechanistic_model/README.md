Required Software
--------

- Python v3.10.12
- Numpy v1.23.5
- SciPy v1.11.4

Usage
--------

For users who want to apply the model to a specific pathway of interest have to create a new class with the same format as `castillom_wsp`, `castillom_sia`, etc. The necessary components are listed below. 

Necessary attributes:

.. code-block:: linux
         
         - self.name # name of the pathway
         - self.index_dgc # index in the y-vector (in the ODE method) that corresponds to DGC (or any other enzyme of interest)
         - self.num_rxns # number of reactions
         - self.num_reactants # number of reactants/species/nodes in the reaction network
         - self.geem_matrix # as explained in code, this is a (num_rxns x num_rxns) matrix, such that g_ij = 1 iff reaction i and reaction j                                  # cannot change simultaneously because they would involve mutations in two separate genetic components. Otherwise                                 # g_ij = 0
         
            
         
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

