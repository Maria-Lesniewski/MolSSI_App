/*
 *  Author: Michael DeLyser
 *
 *  Currently, this computes the rdf for:
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

#include "safe_mem.h"
#include "io_read.h"
#include "gromacs_topology.h"
#include "LDD.h"
#include "rdf_func.h"

#define min2(a,b) (a < b ? a : b)
#define min3(a,b,c) (min2(min2(a,b),c)) 

int main(int argc, char * argv[])
{
  int np, local_rank; // Initiallized some variables to help us keep track of which process we're in in the parallel code
  MPI_Init (&argc, &argv);// From here on out everything is done in parallel [each line is executed by an individual processor]
  MPI_Comm_rank (MPI_COMM_WORLD, &local_rank); // Finds out which process we are in and stores it in local rank
  MPI_Comm_size (MPI_COMM_WORLD, &np); // Finds out how many processes there are total and stores it in np (MPI_COMM_WORLD is a constant representing a predefined communicator that represents all processes)

  int i, j, k;
  double r;
  tW_gmx_trxframe *fr;
  tW_gmx_topology *top;
  top = init_tW_gmx_topology();
  fr = init_tW_gmx_trxframe();
  tW_word top_fnm, fr_fnm;
  tW_RDF_selection *sel;

  int n_params = 10;
  int n_arg_found;
  t_pargs *params = (t_pargs *) ecalloc(n_params,sizeof(t_pargs));
  
  init_arg(&(params[0]),"-s",etSTRING,"input topology filename (req.)",TRUE);
  init_arg(&(params[1]),"-f",etSTRING,"Tinput trajectory filename (req.)",TRUE);
  init_arg(&(params[2]),"-n",etSTRING,"input index filename (opt.)",FALSE);
  init_arg_def(&(params[3]),"-bw",etREAL,"bin width (nm)","0.002");
  init_arg_def(&(params[4]),"-nfr",etINT,"Number of frames (-1 for all)","-1");
  init_arg(&(params[5]),"-NVT",etBOOL,"Signifies NVT simulation",FALSE);
  init_arg(&(params[6]),"-sep",etBOOL,"Save each frame's rdf in its own file",FALSE);
  init_arg(&(params[7]),"-allp",etBOOL,"Compute RDFs between all pairs of atom types",FALSE);
  init_arg(&(params[8]),"-nrdf",etINT,"Number of RDFs to compute",FALSE);
  init_arg(&(params[9]),"-delH",etSTRING,"Input energy fluctutations to compute RDF Derivatives (opt.)",FALSE);

  n_arg_found = get_command_line_args(argc, argv, n_params, params);

  if (local_rank == 0) { print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)); }
  if (n_arg_found <= 0 )
  {
    MPI_Finalize();
    return 0;
  }
  double bw = atof(params[3].value);
  int n_max_frames = atoi(params[4].value);
  if (n_max_frames == -1) { n_max_frames = 1000000000; }
  bool scale_now = (! params[5].bSet); 
  bool save_ind_files = (params[6].bSet);
  bool do_all_pairs = (params[7].bSet);
  int n_rdfs;
  bool derv = (params[9].bSet);

  // scale every frame if you save each frame's RDF
  if (save_ind_files) { scale_now = TRUE; }

  /* Check for mandatory args */
  check_mand_args(params,n_params);

  if ((! params[7].bSet) && (! params[8].bSet))
  {
    fprintf(stderr,"ERROR: either tell me to calculate RDFs between all pairs of atom types with -allp\n");
    fprintf(stderr,"\tOR tell me how many RDFs to calculate with -nrdf\n");
    return 1;
  }

  if ((params[7].bSet) && (params[8].bSet))
  {
    fprintf(stderr,"WARNING: You told me to calculate RDFs between all pairs of atom types\n");
    fprintf(stderr,"\tBut also specified -nrdf. I am ignoring -nrdf\n");
    params[8].bSet = FALSE;
  }

  /* Construct filenames */
  build_filename(top_fnm, params[0].value, eTOP, ".btp");
  build_filename(fr_fnm, params[1].value, eTRAJ, ".trr");

  if (! read_topology(top,(const char *) top_fnm))
  {
    fprintf(stderr,"ERROR: unable to read topology file: %s \n",top_fnm);
    return 1;
  } 

  if (params[2].bSet) { read_index_file(params[2].value,top); }
  generate_generic_index_file(top,params[2].bSet); 

  if (do_all_pairs) 
  { 
    n_rdfs = (top->n_atomtypes * top->n_atomtypes + top->n_atomtypes) / 2; 
  }
  else
  {
    n_rdfs = atoi(params[8].value);
  }
  sel = (tW_RDF_selection *) ecalloc(n_rdfs, sizeof(tW_RDF_selection));
  if (do_all_pairs)
  {
    int rdf_idx = 0;
    for (i = 0; i < top->n_atomtypes; ++i)
    {
      for (j = i; j < top->n_atomtypes; ++j)
      {
        sel[rdf_idx].typeA = i+1+top->contents->mols.nr;
        sel[rdf_idx].typeB = j+1+top->contents->mols.nr;
        ++rdf_idx;
      }
    } 
  }
  else
  {
    int *packed_selections = (int *) ecalloc(n_rdfs * 2, sizeof(int));
    if (local_rank == 0) 
    { 
      get_rdf_selections(top,sel,n_rdfs);
      for (i = 0; i < n_rdfs; ++i)
      {
        packed_selections[2*i] = sel[i].typeA;
        packed_selections[2*i+1] = sel[i].typeB;
      }   
    }
    MPI_Bcast(packed_selections,n_rdfs*2,MPI_INT,0,MPI_COMM_WORLD);// This sends data about the selection to all of our processes, we tell this function that local_rank 0 initialized this data
    if (local_rank != 0) // This initializes the data type for the other processes too
    {
      for (i = 0; i < n_rdfs; ++i)
      {
        sel[i].typeA = packed_selections[2*i];
        sel[i].typeB = packed_selections[2*i+1];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD); // This is a wait order, to make sure that all of the processes get the rdf selection info before we go 
  }
  
  read_first_frame(fr, (const char *) fr_fnm);
  /* These two are necessary because reading the first frame for trr filetype doesn't read the box */
  read_next_frame(fr,FALSE);
  rewind(fr->fp);
  fr->counter = 0;
  /* Now that first frame has been read, we can set up the density profile stuff */
  double long_length = min3(fr->contents->box[0][0],fr->contents->box[1][1],fr->contents->box[2][2])/2.0;
  int n_bins = (int) (long_length / bw) + 1;
  double **rdfs = (double **) ecalloc(n_rdfs, sizeof(double *));
  double **total_rdfs = (double **) ecalloc(n_rdfs, sizeof(double *));
  
  /* Structures required for rdf derivatives */ 
  double **rdfs_d;
  double **total_rdfs_d;
  double *e_vector;

  /* Allocate memory for RDF derivatives calculated from W. Thompson paper, J Chem Phys. 2020 - RJS 01/29/2020 */
  if (derv)
  { 
    rdfs_d = (double **) ecalloc(n_rdfs, sizeof(double *));
    total_rdfs_d = (double **) ecalloc(n_rdfs, sizeof(double *));
  } 

  int n_local_frames = 0, n_total_frames = 0;
  for (i = 0; i < n_rdfs; ++i)
  {
    rdfs[i] = (double *) ecalloc(n_bins, sizeof(double));
    total_rdfs[i] = (double *) ecalloc(n_bins, sizeof(double));

    if (derv)
    {
      rdfs_d[i] = (double *) ecalloc(n_bins, sizeof(double));
      total_rdfs_d[i] = (double *) ecalloc(n_bins, sizeof(double));
    }
  }
  
  int itype;

/* Read in energy fluctuation file, copied code from EM - RJS 01/29/2020 */  
if (derv)
  {
   e_vector = (double *) ecalloc(fr->n_frames, sizeof(double));

    FILE *evec_ptr = fopen(params[9].value, "r");

    int test_sscanf;
    float idx, fr_potl;
    tW_line inp_line;
    int curr_fr = 0;

    while (curr_fr < fr->n_frames)
    {
      get_next_line(evec_ptr,inp_line);
      if ((inp_line[0] != '#') && (inp_line[0] != '@') && (inp_line[0] != '&'))
      {
        test_sscanf = sscanf(inp_line, " %f %f ", &idx, &fr_potl);
        if (test_sscanf != 2)
        {
          fprintf(stderr,"ERROR: unable to read time and energy fluctuations from energy file\n");
          fprintf(stderr,"input line: %i", (int) idx);
          exit(1);
        }
        e_vector[curr_fr] = fr_potl;
        curr_fr++;
      }
    }
  
    fclose(evec_ptr);

  FILE *fpss = open_file("checkner.dat", 'w');
  for (i = 0; i < fr->n_frames; i++)
  {
    fprintf(fpss, "%g \n", e_vector[i]);
  }
  fclose(fpss);
  }

  while ((read_next_frame(fr,FALSE)) && (n_total_frames < n_max_frames))
  {
    ++n_total_frames;
   printf("Popuating rdfs for frame: %d \r", n_total_frames);
    if (n_total_frames % np == local_rank) // This line breaks up the frames read by each processor (e.g. 10 frames, 2 processors, evens will be on 0, odds on 1) 
    {
      ++n_local_frames; // This keeps count of the frames just the one processor has read 
      double ener = 0;
      if (derv)
        ener = e_vector[n_total_frames - 1];
      do_rdfs(top,fr,rdfs,rdfs_d,sel,n_rdfs,bw,n_bins,scale_now, derv, ener);
      if (save_ind_files) 
      { 
        save_rdfs(rdfs, top, sel, n_rdfs, n_bins, bw, n_total_frames, FALSE); 
        save_rdfs(rdfs_d, top, sel, n_rdfs, n_bins, bw, n_total_frames, FALSE);
      }
      if (scale_now) { increment_total_rdfs(total_rdfs, rdfs, total_rdfs_d, rdfs_d, n_local_frames, n_rdfs, n_bins, derv); }
    }
  }
  fclose(fr->fp);
  if ( scale_now )
  {
    for (i = 0; i < n_rdfs; ++i)
    {
      for (j = 0; j < n_bins; ++j)
      {
        rdfs[i][j] = total_rdfs[i][j] * ((double)(n_local_frames) / (double)(n_total_frames)); //We need to cut the weight of each processor on the data
        if (derv) { rdfs_d[i][j] = total_rdfs_d[i][j] * ((double)(n_local_frames) / (double)(n_total_frames)); }
      }
    }
  }
  else
  {
    for (i = 0; i < n_rdfs; ++i)
    {
      for (j = 0; j < n_bins; ++j)
      {
        rdfs[i][j] /= (double)(n_total_frames); // Same as prior before we put everything together we need to weight it by the fraction of frames each processor handled 
        if (derv) { rdfs_d[i][j] /= (double)(n_total_frames); }
      }
    }
  }
  for (i = 0; i < n_rdfs; ++i)
  {
    MPI_Reduce(rdfs[i], total_rdfs[i], n_bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//This says we'd like to handle the rdf data structures by 
    if (derv) { MPI_Reduce(rdfs_d[i], total_rdfs_d[i], n_bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); }
  }
   
  if (local_rank == 0)
  {
    if (! scale_now)
    { 
      scale_rdfs(top, fr, total_rdfs, total_rdfs_d, sel, n_rdfs, bw, n_bins, derv);

/*      double volume = fr->contents->box[0][0] * fr->contents->box[1][1] * fr->contents->box[2][2];
      double factor, dV;
      int rdf_idx;
      for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
      {
        i = sel[rdf_idx].typeA;
        j = sel[rdf_idx].typeB;
        if (i == j) 
        { 
          int hold = top->idx_file->n_at_per_group[i];
          factor = volume / ((double)(hold * (hold) / 2)); 
        }
        else 
        { 
          factor = volume / ((double) (top->idx_file->n_at_per_group[i] * top->idx_file->n_at_per_group[j])); 
        }
        for ( k = 0; k < n_bins; ++k)
        {
          r = ((double)k * bw);
          dV = ((double) (4.0 * M_PI * r * r * bw));
          total_rdfs[rdf_idx][k] *= factor / dV;
        }
      } */
    }
    save_rdfs(total_rdfs, top, sel, n_rdfs, n_bins, bw, n_total_frames, TRUE);
    if (derv) { save_rdfs_d(total_rdfs_d, top, sel, n_rdfs, n_bins, bw, n_total_frames, TRUE); }
  } 

  MPI_Finalize();
  return 0;
}
