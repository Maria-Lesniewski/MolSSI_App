cmake_minimum_required( VERSION 3.0 )
project( "BOCS-FM" )

set ( BOCS_FM_VERSION_MAJOR 2 )
set ( BOCS_FM_VERSION_MINOR 0 )

if( COMMAND cmake_policy )
  cmake_policy(SET CMP0003 NEW)
endif( COMMAND cmake_policy )

set ( CMAKE_BUILD_TYPE Release )

# Change default output directories so not everything lives in src
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# Set compiler flags
set ( CMAKE_C_FLAGS "-O2 -g -fopenmp" )
set ( CMAKE_C_FLAGS_DEBUG "-g" )

# Uncomment this line to use MKL
set ( USE_MKL "ON" )
# And then make sure the path set in line 29 corresponds to the location of your MKL library

set ( RPC_HOME ${RPC_HOME} /data/apps/extern/libtirpc/1.3.2/ ) #Set the location of the parent tirpc (Has lib/share/easybuild stuff) # CHANGE THIS FOR YOUR CLUSTER
include_directories(${RPC_HOME}/include/tirpc) #link to the folder that has rpc and rpcsvc 
set ( RPC_LIB ${RPC_LIB} ${RPC_HOME}/lib/libtirpc.so )

##############################################################################################
# Try to find specified linear algebra libraries 
##############################################################################################
if (USE_MKL STREQUAL "ON")
  if (NOT MKL_LIB)
    set ( MKL_LIB ${MKL_LIB} /data/apps/linux-centos8-cascadelake/gcc-9.3.0/intel-mkl-2020.3.279-omwixzxqb2mkv4hwj3s4pxztqaxzasct/mkl/lib/intel64 )
    # This is the MKL library I tend  to use on the JHU Rockfish cluster, change the path above to your mkl lib [must contain libmkl_rt.so] # CHANGE THIS FOR YOUR CLUSTER
  endif(NOT MKL_LIB)

  find_library( MKLRT mkl_rt ${MKL_LIB} )
  find_library( PTHREAD pthread ${MKL_LIB} )

  if(NOT MKLRT)
    message(FATAL_ERROR "Intel MKL library libmkl_rt.so not found. Check that this file is present at the location specified in your argument to -DMKL_LIB, or try compiling with LAPACK/BLAS instead.")
  endif(NOT MKLRT)

  if(NOT PTHREAD)
    message(FATAL_ERROR "pthreads library libpthreads.so not found. Is it installed? You can point to a nonstandard install location with -DMKL_LIB in your CMake command.")
  endif(NOT PTHREAD)

  set ( LINALG_LIBS ${LINALG_LIBS} ${MKLRT} ${PTHREAD} )

else (USE_MKL STREQUAL "ON") 
  if (NOT LAPACK_LIB)
    set ( LAPACK_LIB ${LAPACK_LIB} /usr/lib64 )
  endif (NOT LAPACK_LIB)
  if (NOT BLAS_LIB)
    set ( BLAS_LIB ${BLAS_LIB} /usr/lib64/ )
  endif (NOT BLAS_LIB)

  find_library(LAPACK lapack ${LAPACK_LIB})
  find_library( PTHREAD pthread ${LAPACK_LIB} )
  find_library(BLAS blas ${BLAS_LIB})

  if(NOT LAPACK)
    message(FATAL_ERROR "LAPACK library liblapack.so not found. Is it installed? You can point to a nonstandard install location with -DLAPACK_LIB in your CMake command.")
  endif(NOT LAPACK)

  if(NOT PTHREAD)
    message(FATAL_ERROR "pthread library libpthread.so not found. Is it installed? You can point to a nonstandard install location with -DLAPACK_LIB in your CMake command.")
  endif(NOT PTHREAD)

  if(NOT BLAS)
    message(FATAL_ERROR "BLAS library libblas.so not found. Is it installed? You can point to a nonstandard install location with -DBLAS_LIB in your CMake command.")
  endif(NOT BLAS)

  set ( LINALG_LIBS ${LINALG_LIBS} ${LAPACK} ${BLAS} )

endif (USE_MKL STREQUAL "ON")


######################################################################################################################

# Set the library combinations we need for compiling different executables
set ( TLIBS ${TLIBS} ${RPC_LIB} m dl )
set ( MLIBS ${MLIBS} ${GLIBS} ${LINALG_LIBS} ${RPC_LIB} m dl )

add_subdirectory( src )

file( COPY src DESTINATION ${CMAKE_BINARY_DIR} ) 
#file( COPY scripts DESTINATION ${CMAKE_BINARY_DIR} ) 

