# BOCS [Add ons for v5]

Contains cmake files and partial BOCS_v5 src code for 2 functions:

1) make_par, an executable I am working on for BOCS_v5 to generate the par.txt file required for forcematching in BOCS based on CG topology info from bocs .btp file [files included to try]

2) o_checker, an executable I wrote as a one off to determine dot products between atomistic_dioxane orientation as a function of the 1-site CG model pair distribution. (e.g. dot product between O-O vectors in atomistic dioxane molecules as a function of CoM distance)
   [files included to try it]

# BOCS Dependencies 
+ gcc/g++ >= 4.9.2 OR icc/icpc >= 13.1
+ CMake >= 3.0
+ OpenMPI >= 1.9.1
+ A linear algebra library, such as LAPACK/BLAS or Intel MKL

