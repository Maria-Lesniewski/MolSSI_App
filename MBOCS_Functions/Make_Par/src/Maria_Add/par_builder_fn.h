/*
 *@file par_builder_fn.h
 *@authors Maria Lesniewski 
 *@brief Helper functions for the buildpar executable [driver = par_builder.c]
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../BOCS_Default/cgff_types.h"
#include "../BOCS_Default/gromacs_topology.h"
#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/wnoid_math.h"

enum { mNB, mBOND, mANGLE, mDIH, mIMNB}; // used for our wrapper functions to determine which top we're listing through, add to the switches in the wrapper and this as types are added
void write_par(tW_word fnm, tW_gmx_topology *top); // main function associated with this header

// Functions for printing interaction lists in par.txt format
void print_imnb( FILE *fp, tW_word** ilist, int len);
void print_dihedral( FILE *fp, tW_word** ilist, int len);
void print_angle( FILE *fp, tW_word** ilist, int len);
void print_bond( FILE *fp, tW_word** ilist, int len);
void print_pair( FILE *fp, tW_word** ilist, int len);
void write_interactions(FILE *fp, tW_word** ilist, int len, int mode);
void total_print_imnb( FILE *fp, tW_word** ilist, int len);
void total_print_dihedral( FILE *fp, tW_word** ilist, int len);
void total_print_angle( FILE *fp, tW_word** ilist, int len);
void total_print_bond( FILE *fp, tW_word** ilist, int len);
void total_print_pair( FILE *fp, tW_word** ilist, int len);
void write_total_list( FILE *fp, tW_word** ilist, int len, int mode);

// Functions for listing required interaction types while looping through system top
void fill_impropernb(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check);
void fill_impropernb_list(tW_gmx_topology* top, tW_word **ilist, int list_len);
void fill_dihedrals(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check);
void fill_dihedral_list(tW_gmx_topology* top, tW_word **ilist, int list_len);
void fill_angles(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check);
void fill_angle_list(tW_gmx_topology* top, tW_word **ilist, int list_len); 
void fill_bonds(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check); 
void fill_bond_list(tW_gmx_topology* top, tW_word **ilist, int list_len);  
int fill_intermol_pairs(tW_word **ilist, tW_gmx_topology* top, int ifoundsofar);
int molecule_choose_2_interactions(tW_word **ilist, tW_molecule *mol, int ifoundsofar);
int fill_nb_list_unique_mol(tW_word **ilist, tW_molecule *mol, int ifoundsofar);
void brutal_list_pair_int(tW_gmx_topology* top, tW_word **ilist, int list_len);
void list_choose_2_interactions(tW_gmx_topology* top, tW_word **ilist, int list_len);
void fill_nb_list(tW_gmx_topology* top, tW_word **ilist, int list_len);

//Functions for counting required interaction types while looping through system top
int search_improper(tW_molecule * mol, bool * beq_found, int check);
int get_n_impropernb(tW_gmx_topology* top);
int list_intermol_interactions(tW_word ** ilist, tW_gmx_topology* top, int* len_ilist);
int list_cross_interactions(tW_word **ilist, tW_molecule* mol, int* len_ilist);
int list_self_interactions(tW_word **ilist, tW_molecule* mol, int* len_ilist);
bool check_pairlist(tW_word typeA, tW_word typeB, tW_word **ilist, int len_list);
bool check_list(int* list_of_int, int len, int my_int);
int add_to_nb_list(tW_word ** ilist, tW_molecule * mol, int* len_ilist);
void sort_moltypes(tW_gmx_topology *top, int* one_list, int* multi_list, int* dim_holder);
int brutal_count_nb_interactions(tW_gmx_topology *top);
bool naive_exclusion_check(tW_gmx_topology *top);
int count_required_nb(tW_gmx_topology *top);
int search_dihedrals(tW_molecule * mol, bool * beq_found, int check);
int get_n_distinct_dihedrals(tW_gmx_topology *top);
int search_angles(tW_molecule * mol, bool * beq_found, int check);
int get_n_distinct_angles(tW_gmx_topology *top);
int get_n_distinct_bonds(tW_gmx_topology *top); // fill_bond_list, but count the unique bonds
int search_bonds(tW_molecule * mol, bool * beq_found, int check); // fill_bonds, but count them
int choose(int n, int k);
int calc_total_interactions(tW_gmx_topology *top);
