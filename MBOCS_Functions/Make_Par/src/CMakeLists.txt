# 1) Automate finding BLAS/LAPACK or MKL with provided modules
# 2) Can I automate a find gromacs, or have the path location to gmxdir
# available as a command line argument?




add_executable (make_par Maria_Add/par_builder.c
			 Maria_Add/par_builder_fn.c
			 BOCS_Default/safe_mem.c
			 BOCS_Default/wnoid_math.c
			 BOCS_Default/io_read.c
			 BOCS_Default/gromacs_topology.c )

target_link_libraries( make_par ${TLIBS} )
