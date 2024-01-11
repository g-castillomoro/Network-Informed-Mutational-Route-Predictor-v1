Required Software
--------

- Python v3.10.12
- Numpy v1.23.5
- SciPy v1.11.4

Usage
--------

For users who want to apply the model to a specific pathway of interest have to create a new class with the same attributes as `castillom_wsp`, `castillom_sia`, etc. The necessary components are listed below. 

Necessary attributes:
.. code-block:: linux
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

