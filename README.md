# FenicsFEA
Describe the cooling time and the temperature field for a given material slab.

2d_solver: main solver call

light_abs_DT : calculates the temperature increase for the slab when light of a given wavelenght is shone on it. It wants a .csv file with the n&k values,depending on the material modelled

max_q_outflow_bulk : calculates the max outflow for the current slab temperature. This is the Neuman BC for the weak formulation of the problem.

