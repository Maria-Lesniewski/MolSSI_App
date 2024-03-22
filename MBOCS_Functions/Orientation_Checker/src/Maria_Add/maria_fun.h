#ifndef MARIA_FUN_H
#define MARIA_FUN_H

#ifdef __cplusplus 
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"

#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/gromacs_topology.h"

//This is just a more generally described rdf selection structure from rdf.c by MRD&RJS
typedef struct tW_user_selection{
	int sel_one;
	int sel_two;
} tW_user_selection;

typedef struct tW_molcom_vector{
 	double com_position[3]; // This will hold the molecule position after mapping to 1 site
	double orientation_vector[3]; // This will hold the orientation vector specified 
} tW_molcom_vector;


//These functions are altered versions of the rdf_func.c package
double get_rect_prism_dx(dvec xi, dvec xj, tW_matrix box, dvec dr); // not this one though, we ripped this one from rdf_fun.c to avoid an extra include
void get_user_selections(tW_gmx_topology *top, tW_user_selection *sel);
bool check_valid_index_selections(tW_index_file *idx_file, int c1, int c2);
int get_idx(double dist, double bw, int n_bins);

//These functions are custom to this script, more detailed doc is in the accompanying .c file
double get_mol_vect(dvec xi, dvec xj, tW_matrix box, dvec dr);
void calc_molcom_vector(tW_gmx_trxframe *fr, tW_gmx_trxframe *cgfr, tW_gmx_topology *top, tW_molcom_vector* comv, int* packed_selections);
void count_dot_dist(tW_molcom_vector* comv, double* dot_dis_storage, int N_sites, tW_matrix box, double bw, int n_bins);
void store_molcom_from_cgframe(tW_gmx_trxframe *cgfr, tW_molcom_vector* comv);
void store_vector_from_aatop(tW_gmx_trxframe *fr, tW_gmx_topology *top, tW_molcom_vector* comv, int* packed_selections);
double lp2_dot(double* vect_i, double* vect_j);
double legendre_p_second_rank(double x);
double get_dot_product(dvec vect_i, dvec vect_j);
void write_binned_vector(double *vector, int n_bins, double bw, char* filename);
double simple_scale(double number_to_scale, double bw, int idx);
double density_scale(double number_to_scale, tW_matrix box, int N_sites);
void alt_count_dot_dist(tW_molcom_vector* comv, double* dot_dis_storage, int N_sites, tW_matrix box, double bw, int n_bins, int * bin_counter);
#ifdef __cplusplus
}
#endif
#endif
