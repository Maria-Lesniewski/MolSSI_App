/*
 *  Author: Maria Lesniewski
 *  working off old RDF script by R. Szukalo and M. Delyser
 * Documentation to change
 * Currently, this computes the rdf for:
 *  -- all pairs of atomtypes if -allp specified
 *  -- only specified pairs if -nrdf # is specified
 *  There are no bonded exclusions
 *  It works great for CG trajectories, though
 *  It is parallelized with MPI
 *  
 *  Note, parallel works on aci-b for batch jobs,
 *  but is sometimes crashes when running in parallel
 *  on interactive jobs (some of the processors segfault
 *  when trying to read the first number from the .trr file)
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "../BOCS_Default/read_map.h"
#include "../BOCS_Default/transf_map.h"
#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/gromacs_topology.h"
#include "maria_fun.h"

#define min2(a,b) (a < b ? a : b) // ML returns the min value of 2 values passed
#define min3(a,b,c) (min2(min2(a,b),c)) // The same but with three

int main(int argc, char * argv[])
{
//This stuff will be used for MPI Processing
  int np, local_rank;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &local_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np);

//Initializing BOCS data structures 
  FILE *fp_map; 
  tW_site_map *CG_map;
  tW_gmx_trxframe *fr, *cgfr; //The traj we're reading and one we're making for CoM Mapping 
  tW_gmx_topology *top, *cgtop; // The aa btp we're gonna read from and the cg_file we'll write to 

  top = init_tW_gmx_topology(); 
  cgtop = init_tW_gmx_topology();
  fr = init_tW_gmx_trxframe();
  cgfr = init_tW_gmx_trxframe();
  
  tW_word top_fnm, fr_fnm; // This is so that we can store the names of our structures
  tW_user_selection *sel; // This is what I would call a tworay that stores ints representing the selection

//Initilizing other Variables we need
  int N_sites; 

//Initializing Space for User Input / Error Handling / Take Args from Command Line
  int n_params = 6;
  int n_arg_found;
  t_pargs *params = (t_pargs *) ecalloc(n_params,sizeof(t_pargs)); // Leave room for the user's input 6 things needed
  
  init_arg(&(params[0]),"-s",etSTRING,"input topology filename (req.)",TRUE); // 1) Pass the top file with a -s, this just creates that slot 
  init_arg(&(params[1]),"-f",etSTRING,"input trajectory filename (req.)",TRUE); // 2) Pass the trajectory with a -f
  init_arg(&(params[2]),"-n",etSTRING,"input index filename (opt.)",FALSE); // 3) Pass the index with a -n
  init_arg(&(params[3]),"-p",etSTRING,"input mapping topology file (req.)",TRUE); //4)Pass the mapping topology file -p
  init_arg_def(&(params[4]),"-bw",etREAL,"bin width (nm)","0.002"); // 5) This sets the binsize we're gonna use 
  init_arg_def(&(params[5]),"-nfr",etINT,"Number of frames (-1 for all)","-1"); //  6) This sets how many frames we're gonna read

  n_arg_found = get_command_line_args(argc, argv, n_params, params); // Shove user lines into data structures for file names 

//We only want the root processor to print this table while troubleshooting user input 
  if (local_rank == 0) { print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)); }
  if (n_arg_found <= 0 )
  {
    MPI_Finalize();
    return 0;
  }

//Convert User Input to appropriate data type
  double bw = atof(params[4].value);
  int n_max_frames = atoi(params[5].value);
  if (n_max_frames == -1) { n_max_frames = 1000000000; }

  /* Check for mandatory args */
  check_mand_args(params,n_params);


  /* Construct filenames */
  build_filename(top_fnm, params[0].value, eTOP, ".btp");// This fills the top_fnm string with what file name we typed if it passes some checks 
  build_filename(fr_fnm, params[1].value, eTRAJ, ".trr"); // Likewise for fr_fnm 

  if (! read_topology(top,(const char *) top_fnm)) // Sets the kind of topology we're using based on extension name, calls the appropriate reader in to fill the tW_top structure, returns false if type not detected
  {
    fprintf(stderr,"ERROR: unable to read topology file: %s \n",top_fnm);
    return 1;
  } 

  if (params[2].bSet) { read_index_file(params[2].value,top); }
  generate_generic_index_file(top,params[2].bSet); 
// Finish reading files off command line, now actively ask user for info 

  sel = (tW_user_selection *) ecalloc(1, sizeof(tW_user_selection));
  int *packed_selections = (int *) ecalloc(2, sizeof(int));
    if (local_rank == 0) 
    { 
      get_user_selections(top,sel);
      packed_selections[0] = sel->sel_one;
      packed_selections[1] = sel->sel_two;    
    }
    MPI_Bcast(packed_selections,2,MPI_INT,0,MPI_COMM_WORLD);
    if (local_rank != 0)
    {
        sel->sel_one = packed_selections[0];
        sel->sel_two = packed_selections[1];
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (local_rank == 0) {
    printf("Congratulations, we've taken user input and packed it\n");
    }

 if (top->idx_file->n_at_per_group[packed_selections[0]] != top->idx_file->n_at_per_group[packed_selections[1]]){
 printf("You Can't Use These Selections For Vectors because there aren't the same qty\n");
 printf("There are %d atoms in sel1 and %d atoms in sel2", top->idx_file->n_at_per_group[packed_selections[0]], top->idx_file->n_at_per_group[packed_selections[1]]);
 exit(1);
 }

// Now Onto Analysis Specific Set Up
  read_first_frame(fr, (const char *) fr_fnm);
  printf("we read the first frame\n");
  /* These two are necessary because reading the first frame for trr filetype doesn't read the box */
  read_next_frame(fr,TRUE);
  printf("we read the next frame\n");
  rewind(fr->fp);
  printf("we rewinded\n");
  fr->counter = 0;

/*read CG mapping */
  fp_map = open_file(params[3].value, 'r');
  CG_map = get_CG_map(fp_map, &N_sites);
  printf("Read map\n");
/* Use mapping and AA frame to build CG frame structure */
  
  setup_CG_trr(N_sites, cgtop, fr, cgfr);
  printf("set up the CG trr struct\n");
  printf("Lets try to set up an index file for the top\n");
  generate_generic_index_file(cgtop, FALSE);
  printf("generated index file\n");
/* Map the first atomistic frame */ 
 map_config(N_sites, CG_map, fr, cgfr);
  printf("mapped the first frame\n");
// Now that the first frame has been read, we can set up the required box info for our pbc and frame dependent structures
  double pbc_cutoff_length = min3(fr->contents->box[0][0],fr->contents->box[1][1],fr->contents->box[2][2])/2.0;
  printf("box lengths: %f %f %f \n", fr->contents->box[0][0], fr->contents->box[1][1],fr->contents->box[2][2]);

  int n_bins = (int) (pbc_cutoff_length / bw) + 1;
  printf("We're using %i bins based on half of the short box length, the total length will be %f of the first frame, user beware if you are passing a non-equillibrated NPT simulation\n",\
		  n_bins, (double)(n_bins*bw));

  double* dot_dis = (double *) ecalloc(n_bins, sizeof(double)); // This is going to hold our results from adding dot products by CoM pair distances 
  double* t_dot_dis = (double *) ecalloc(n_bins, sizeof(double)); // This is going to hold the merged results from different threads
  int* bin_counter = (int *) ecalloc(n_bins, sizeof(int)); // This is to store a count of things going in the bin if we use alt_dot_dis in the analyis
  int* t_bin_counter = (int *) ecalloc(n_bins, sizeof(int)); // This will hold the total counter post thread merge
  printf("Allocated Memory for dot dis\n");
 tW_molcom_vector* com_orientation_frame_info = (tW_molcom_vector *) ecalloc(N_sites, sizeof(tW_molcom_vector));
 printf("Allocated Memory for %d different molecules per frame \n", N_sites);
 int n_total_frames = 0;
 int n_local_frames = 0;

 // Let's Do the Actual Analysis Now
 
 while ((read_next_frame(fr,FALSE)) && (n_total_frames < n_max_frames)){
	 ++n_total_frames;
	 printf("Analyzing frame %d \n", n_total_frames);
	if (n_total_frames % np == local_rank){ // This line breaks up the frames read by each processor (e.g. 10 frames, 2 processors, evens will be on 0, odds on 1) 
		++n_local_frames;
		map_config(N_sites, CG_map, fr, cgfr);
		calc_molcom_vector(fr, cgfr, top, com_orientation_frame_info, packed_selections);
		//count_dot_dist(com_orientation_frame_info, dot_dis, N_sites, fr->contents->box, bw, n_bins);
		alt_count_dot_dist(com_orientation_frame_info, dot_dis, N_sites, fr->contents->box, bw, n_bins, bin_counter);
	}	
 }
 fclose(fr->fp);
// Now we need to scale everything down so that we can merge processes with weight according to the qty of data handled
for (int i = 0; i < n_bins; i++){
	//if (i==400) {printf("MERGING bin %d with total count %d \n", i, bin_counter[i]);
	//printf("Value before merge: %f", dot_dis[i]);}
	//dot_dis[i]= dot_dis[i] * ((double)(n_local_frames)/((double)(n_total_frames)* N_sites));
	dot_dis[i] = dot_dis[i] * ((double)(n_local_frames)/((double)(n_total_frames)));
	//if (i==400){printf("Value after: %f should be 865 times smaller and max_frames %f and total_frames %f\n and total_bins should be same : %d \n", dot_dis[i], n_max_frames, n_total_frames, bin_counter[i]);}
}
printf("bin_count before merge: %d \n", bin_counter[400]);
// MERGE DATA
MPI_Reduce(dot_dis, t_dot_dis, n_bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Reduce(bin_counter, t_bin_counter, n_bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
printf("bin 400 post merge %d \n", t_bin_counter[400]);

//Average
for (int j = 0; j < n_bins; j++){
	//t_dot_dis[j] = t_dot_dis[j]/n_total_frames; // Use this line when using count_dot_dis in the Actual Analysis Loop
	if (t_bin_counter[j] != 0){
	t_dot_dis[j] = t_dot_dis[j] / (double)(t_bin_counter[j]);}
	if (j ==400){
	printf("Hey jus in bin 400 again, before we divide \n");}
}
//Save Results, we only need a copy from the main thread
if (local_rank == 0 ){
write_binned_vector(t_dot_dis, n_bins, bw, "Dot_Distance_Data.dat");
}
 printf("Processor %d made it to the end somehow\n", local_rank);
MPI_Finalize();
return 0;
}
