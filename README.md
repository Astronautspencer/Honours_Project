
# Honours_Project


# The role of mesh when fitting LGCP using INLA-SPDE approach - simulation study

This repository contains the results and source code for the simulation study for the Honours project.


#### Main files

- `source_code/help_functions.R` Functions containing the implementation of all the differnet approximation methods for the integral in question, in addition to a series of help functions related to the simulation experiments.
- `source_code/full_simulations_script.R` The script used to carry out the full simulation experiment in the paper. The script also contains the short test version that quickly runs trough some parameter combinations (to check that the script works as intended).
- `source_code/permutation_tests.R` The script used to carry out the permutation tests mentioned in the paper.
- `source_code/illustrating_mesh_approximated_field.R` A self-contained script for illustating how to approximate a given field on specific mesh using the finite element method (FEM).





#### Full result table for simulations experiments
A sort and searchable table with the full simulation results is available [here](https://martinju.github.io/LGCP-normConst-simulations/sim_res.html).

#### Permutation tests
A sort and searchable table with results for permutation tests of pairwise differences between  all combinations 
of the approximation methods, for each parameter combination is available 
[here](https://martinju.github.io/LGCP-normConst-simulations/permut_tests.html).

