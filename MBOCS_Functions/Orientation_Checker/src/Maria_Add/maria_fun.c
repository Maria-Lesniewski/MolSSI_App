#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"

#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/wnoid_math.h"
#include "../BOCS_Default/gromacs_topology.h"
#include "maria_fun.h"

#define mymin(a,b) (a < b ? a : b)
#define mymax(a,b) (a > b ? a : b)


/*get_rect_prism_dx(position_i, poosition_j, boxfromtraj, absolute displacement vector storage)
 *This is a function originally written by MRD and RJS for the rdf function
 *It receives pointers we get from the trajectory frame and populates an absolute displacement vector
 *This absolute displacement vector accounts for the period images that may be implied by the coordinates
 *The function itself returns the magnitude of this vector, which makes it great for binning by dist. 
 *One note: This function does not maintain the directionality of the displacement under PBC
 */
double get_rect_prism_dx(dvec xi, dvec xj, tW_matrix box, dvec dr)
{
  int d;
  double delta, dr2 = 0.0;
  for (d = 0; d < 3; ++d)
  {
    delta = fabs(xi[d] - xj[d]);
    dr[d] = mymin(delta,box[d][d] - delta);
    //printf("inside the function we have %d i = %f and j = %f \n",d, xi[d], xj[d] );
    dr2 += dr[d] * dr[d];
  }
  //printf("disp before return %f \n", sqrt(dr2));
  return (sqrt(dr2));
}

/*get_mole_vect()
 *This function takes two position vectors (xi xj) and returns the distance (with pbc) 
 *displacement vectors are stored in the new dr vec
 *The dr vector maintains the orientational directionality imposed by subtracting coordinates
 *This function handles pbc with the assumption that the displacement we seek is the intramolecular direction vector
 */
double get_mol_vect(dvec xi, dvec xj, tW_matrix box, dvec dr)
{ // This is tricky because the molecule might be split across pbc and if it is the vector sign will be flipped 
  int d;
  bool broken; 
  double delta, dr2 = 0.0;
  for (d = 0; d < 3; ++d)
  {
    delta = xi[d] - xj[d];
    // check if split by pbc
    if (fabs(delta) > box[d][d] - fabs(delta)){
    printf("warning molecules split across pbc, doing some math to handle it\n");
	dr[d]= -1* (box[d][d] - fabs(delta)); // If the image is split, the displacement vector is backwards
    }
    else{ // If its not split then we'll leave delta as is
    dr[d] = delta; 
    }
    dr2 += dr[d] * dr[d];
  }
  return (sqrt(dr2));
}
//get_user_selections(top: topology structure, sel: blank spot for selection storage)
//This is a simplified version of the function get_RDF_selections written by MRD and RJS
//this function will print options from the passed (from ndx) or infered index (from btp) options. It will then take the options specified and store them
//inside of the sel structure passed 
void get_user_selections(tW_gmx_topology *top,tW_user_selection *sel)
{
  int i, idx1, idx2;
  print_index_names_and_atom_counts(top->idx_file);
    do
    {
      fprintf(stdout, "We need to define an orienation vector based on the positions of atoms in each molecule\n");
      fprintf(stdout, "Please enter two atom types that you'd like a displacement vector between\n");
      fprintf(stdout,"Enter index for selection one  %d: \n",i);
      scanf("%d",&idx1);
      fprintf(stdout,"Enter index for selection two %d: \n",i);
      scanf("%d",&idx2);
    } while ( ! check_valid_index_selections(top->idx_file, idx1, idx2) );
      sel->sel_one = idx1;
      sel->sel_two = idx2;
      fprintf(stdout,"The function will use %s and %s as our selections from the index file \n",top->idx_file->group_names[idx1],top->idx_file->group_names[idx2]);
}

/*check_valid_index_selections(): 
 * Returns false if user puts bad indices in on the selection menu
 * used for trouble shooting, again, ripped straight from rdf_fun.c by MRD and RJS
 */
bool check_valid_index_selections(tW_index_file *idx_file, int c1, int c2)
{
  if (c1 >= idx_file->n_groups)
  {
    fprintf(stderr,"ERROR: selection %d is too high.\n",c1);
    return FALSE;
  }
  if (c2 >= idx_file->n_groups)
  {
    fprintf(stderr,"ERROR: selection %d is too high.\n",c2);
    return FALSE;
  }
  return TRUE;
}

/*calc_molcom_vector(aa_frames, cg_frames, aa_top, molcom_storage_vector)
 *This function receives high resolution and coarse resolution frames as input 
 * It uses the corresponding topologies and user selection to get the com and orientation 
 * vector for each molecule
 * This function makes the rather large assumptions that 
 * 	1) You have 1 pair of atoms you're using as your vector per molecule 
 * 	2) The molecule is mapped to 1 site [you'll have to alter this and the executable to take CG selections from the user if you want more than that]
 * 	Note that this change would be non-trivial in the driver code, I suggest mapping in advance and reading in an index file if you're interested in writing 
 * 	a general function that gives the orientation of multiple kinds of CG sites
 */
void calc_molcom_vector(tW_gmx_trxframe *fr, tW_gmx_trxframe *cgfr, tW_gmx_topology *top, tW_molcom_vector* comv, int* packed_selections){

  	store_molcom_from_cgframe(cgfr, comv);
	store_vector_from_aatop(fr, top, comv, packed_selections);

}
/*store_molcom_from_cgframe(cgframes, comv )
 *Stores the CoM coordinates from the mapped trajectory frame
 *Into the com_position field of the tW_molcom_vector data structure 
 */
void store_molcom_from_cgframe(tW_gmx_trxframe *cgfr, tW_molcom_vector* comv){
	//printf("Inside the store_molcom function\n");
	int total_cg_sites = cgfr->contents->natoms; 
	printf("set total sites to %d\n", total_cg_sites);
	dvec x; 
	for (int i = 0; i < total_cg_sites; i++){
		get_pos_of(cgfr, i, x);
                comv[i].com_position[0]=x[0];
                comv[i].com_position[1]=x[1];
                comv[i].com_position[2]=x[2];
                //printf("MOL ID: %d COM POSITION: %f %f %f \n", i, comv[i].com_position[0], comv[i].com_position[1], comv[i].com_position[2]);
	}
}
/*store_vector_from_aatop(aaframe, aatop, storage_CoM_Orientation vector, packed_user_selections)
 *Obtains the orientation vectors defined by the user selection and stores into the CoM_Orientation vector for the molecule
 *The vector convention is the second selection coordinates minus the first one
 *Note this function takes care to retain the relative directionality of the unit vectors (in case 3D quadrants matter in the intended function)
 */
void store_vector_from_aatop(tW_gmx_trxframe *fr, tW_gmx_topology *top, tW_molcom_vector* comv, int* packed_selections){
	//printf("Inside of the vector store function \n");
	int total_atoms_in_either_group= top->idx_file->n_at_per_group[packed_selections[0]]; // earlier we checked that this matched for the vector, see the driver code 
	//printf("Detected %d total vectors \n", total_atoms_in_either_group);
	int sel1_idx;
	int sel2_idx;
	dvec dr;
	double dr_length;
	for (int i = 0; i < total_atoms_in_either_group; i++){
		sel1_idx=top->idx_file->atoms_in_group[packed_selections[0]][i];
		sel2_idx=top->idx_file->atoms_in_group[packed_selections[1]][i];
		
		dr_length = get_mol_vect(fr->contents->x[sel2_idx], fr->contents->x[sel1_idx], fr->contents->box, dr); // selection2-1
		
		comv[i].orientation_vector[0]=dr[0]/dr_length;
		comv[i].orientation_vector[1]=dr[1]/dr_length;
		comv[i].orientation_vector[2]=dr[2]/dr_length;
		//printf("Orientation vector %d : %f %f %f\n", i, comv[i].orientation_vector[0], comv[i].orientation_vector[1], comv[i].orientation_vector[2] );
	}
}
/*count_dot_dist(comv: Center_of_Mass_Orientation_Vector_Storage, dot_dis_storage: the vector to store statistics in, N_sites: total_sites, box: simulation box (from frame)...
 *               ... bw : bin width for dot_dis_storage, n_bins : len(dot_dis_storage)) 
 * This is the meat of the executable written here. This function finds the CoM displacements and Orientation Dot Product for every unique pair of molecules
 * It then bins the Legendre Polynomial of that dot product by the CoM displacement
 */
void count_dot_dist(tW_molcom_vector* comv, double* dot_dis_storage, int N_sites, tW_matrix box, double bw, int n_bins){
	//printf("Inside the Dot Product Binning Function\n");
	dvec dummydr;	
	int idx;
	double displace;
	double o2l_polynomial_dot;
	double scaled_o2l;
	for (int i = 0; i < N_sites ; i++){
		for (int j = i +1; j < N_sites; j++){
		  //Find the CoM dist (Note these typings work because dvec is a 3 entry double vec)
		  displace=get_rect_prism_dx(comv[i].com_position, comv[j].com_position, box, dummydr);
		  if(displace < 0.4){
		  //printf("CHECK THESE INDICES %d %d \n", i,j );
		  //printf("We passed vector %d: %f %f %f and vector %d: %f %f %f \n", i, comv[i].com_position[0], comv[i].com_position[1], comv[i].com_position[2], j, \
		//		  comv[j].com_position[0], comv[j].com_position[1], comv[j].com_position[2]);
		  //printf("displace %f\n",displace);
		  }
		  //ID a bin index
		  idx=get_idx(displace, bw, n_bins);
		  //if(i==0 && j == 1){printf("We've gone ahead and placed this in bin %d given the bw %f \n", idx, bw);}
		  //if(i==654 && j == 655){printf("We've gone ahead and placed this in bin %d given the bw %f \n", idx, bw);}
		  //Find the 2nd Order Legendre Polynomial of the Angle Between Molecules
		  //1/2(3x^2-1), where x = cos(\theta), and \theta is the angle btwn mol vectors 
		  o2l_polynomial_dot = lp2_dot(comv[i].orientation_vector, comv[j].orientation_vector);
		  //if(i==0 && j==1){printf("We got %f for the o2l legendre polynomial dot product\n", o2l_polynomial_dot);}
		  //if(i==654 && j==655){printf("We got %f for the o2l legendre polynomial dot product\n", o2l_polynomial_dot);}
		  //scale
		  //scaled_o2l = density_scale(simple_scale(o2l_polynomial_dot, bw, idx), box, N_sites);
		  //store
		  //dot_dis_storage[idx] = dot_dis_storage[idx] + o2l_polynomial_dot;
		  dot_dis_storage[idx] = dot_dis_storage[idx] + density_scale(simple_scale(2, bw, idx), box, N_sites); 
		  //dot_dis_storage[idx] = dot_dis_storage[idx] + scaled_o2l;
		  //dot_dis_storage[idx]= dot_dis_storage[idx] + density_scale(simple_scale(1, bw, idx),box, N_sites); 
		}
	}
}
/*get_idx(item_to_bin, binwidth, total_bins)
 *This function was written by MRD and RJS in the rdf_func.c packages for rdf.c
 *As it says- this returns the appropriate index for histograming data given some binsize and length */
int get_idx(double dist, double bw, int n_bins)
{
  int idx = (int)((dist+0.5*bw) / bw);
  if (idx >= n_bins) { idx = -1; }
  return idx;
}
/*lp2_dot (storage_vector, vector_to_dot_prod1, vector_to_dot_prod2)
 *Returns the second rank legendre polynomial of the dot product between two vectors*/
double lp2_dot(double* vect_i, double* vect_j){
	double dot_product = get_dot_product(vect_i, vect_j);
	return legendre_p_second_rank(dot_product);
}
/*legendre_p_second_rank(x, as in p(x))
 *Returns the second order legendre polynomial of an arg*/
double legendre_p_second_rank(double x){
	return 0.5*(3*x*x-1);
}
/*get_dot_product(vector_to_dot1, vector_to_dot2)
 *Returns the dot product of two vectors*/
double get_dot_product(dvec vect_i, dvec vect_j){
	return (vect_i[0]*vect_j[0] + vect_i[1]*vect_j[1] + vect_i[2]*vect_j[2]);
}
/*write_binned_vector(y_data, total_bins, bin_increment(binwidth), stringoutputfilename)
 * Writes an x y file for a vector passed in such that the y column is the contents of the vector
 * and the x is the associated domain as inferred by the bin spacing and vector indices
 */
void write_binned_vector(double *vector, int n_bins, double bw, char* filename){
 	double x_column;
	FILE* fp = open_file(filename, 'w');
 	  for (int i = 0; i < n_bins; i++){
	       	x_column = (double) i * bw;
		fprintf(fp, "%f %f \n", x_column, vector[i]);
	  }
  	fclose(fp);	
}
/*simple_scale(number_to_scale, bin width, index)
 * for a function f(x), returns f(x)/(4PIx^2), sort of like an RDF uses 
 */
double simple_scale(double number_to_scale, double bw, int idx){
	double r = (double) idx * bw;
	double dV = r*r*4.0*M_PI*bw;

	return (number_to_scale/dV);
}
/*density_scale(number_to_scale, box, total_N)
 * Returns a number scaled by the number density 
 */
double density_scale(double number_to_scale, tW_matrix box, int N_sites){
	double number_density = (double) N_sites / (box[0][0]*box[1][1]*box[2][2]);	
	return number_to_scale/number_density;
}
/*alt_count_dot_dist
 *see count_dot_dist documentation, this behaves much the same way except for it counts the total of L2P's added to the bin 
 * and does not scale the qty. like an rdf. The intent is that the entries in bin_counter can be used later for averaging 
 */
void alt_count_dot_dist(tW_molcom_vector* comv, double* dot_dis_storage, int N_sites, tW_matrix box, double bw, int n_bins, int * bin_counter){
        //printf("Inside the Dot Product Binning Function\n");
        dvec dummydr;
        int idx;
        double displace;
        double o2l_polynomial_dot;
        double scaled_o2l;
        for (int i = 0; i < N_sites ; i++){
                for (int j = i +1; j < N_sites; j++){
                  //Find the CoM dist (Note these typings work because dvec is a 3 entry double vec)
                  displace=get_rect_prism_dx(comv[i].com_position, comv[j].com_position, box, dummydr);
                  if(displace < 0.4){
                  //printf("CHECK THESE INDICES %d %d \n", i,j );
                  //printf("We passed vector %d: %f %f %f and vector %d: %f %f %f \n", i, comv[i].com_position[0], comv[i].com_position[1], comv[i].com_position[2], j, \
                                  comv[j].com_position[0], comv[j].com_position[1], comv[j].com_position[2]);
                  //printf("displace %f\n",displace);
                  }
                  //ID a bin index
                  idx=get_idx(displace, bw, n_bins);
                  //Find the 2nd Order Legendre Polynomial of the Angle Between Molecules
                  //1/2(3x^2-1), where x = cos(\theta), and \theta is the angle btwn mol vectors 
                  o2l_polynomial_dot = lp2_dot(comv[i].orientation_vector, comv[j].orientation_vector);
                  //store and count how many things are stored
                  dot_dis_storage[idx] = dot_dis_storage[idx] + o2l_polynomial_dot;
		  bin_counter[idx] = bin_counter[idx] + 1; 
                }
        }
}

