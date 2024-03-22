#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"

#include "safe_mem.h" 
#include "io_read.h"
#include "wnoid_math.h"
#include "gromacs_topology.h"
#include "rdf_func.h"

#define mymin(a,b) (a < b ? a : b)
#define mymax(a,b) (a > b ? a : b)
double get_rect_prism_dx(dvec xi, dvec xj, tW_matrix box, dvec dr)
{
  int d;
  double delta, dr2 = 0.0;
  for (d = 0; d < 3; ++d)
  {
    delta = fabs(xi[d] - xj[d]);
    dr[d] = mymin(delta,box[d][d] - delta);
    dr2 += dr[d] * dr[d];
  }
  return (sqrt(dr2));
}

int get_idx(double dist, double bw, int n_bins)
{
  int idx = (int)((dist+0.5*bw) / bw);
  if (idx >= n_bins) { idx = -1; }
  return idx;
}

bool check_valid_rdf_selections(tW_index_file *idx_file, int c1, int c2)
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
//get_rdf_selections(top: topology structure, sel: blank spot for selection storage, n_rdfs: total number of rdfs to calculate)
//If the user doesn't want all of the possible rdfs to be calculated we need to ask them for which one
//this function will print options from the passed (from ndx) or infered index (from btp) options. It will then take the options specified and store them
//inside of the sel structure passed 
void get_rdf_selections(tW_gmx_topology *top,tW_RDF_selection *sel, int n_rdfs)
{
  int i, idx1, idx2;
  print_index_names_and_atom_counts(top->idx_file);
  for (i = 0; i < n_rdfs; ++i)
  {
    do
    {
      fprintf(stdout,"Enter index for group A of RDF %d: \n",i);
      scanf("%d",&idx1);
      fprintf(stdout,"Enter index for group B of RDF %d: \n",i);
      scanf("%d",&idx2);
    } while ( ! check_valid_rdf_selections(top->idx_file, idx1, idx2) );
      sel[i].typeA = idx1;
      sel[i].typeB = idx2;

      fprintf(stdout,"rdf %d will be between %s and %s \n",i,top->idx_file->group_names[idx1],top->idx_file->group_names[idx2]);

  }
}

void scale_rdfs(tW_gmx_topology *top, tW_gmx_trxframe *fr, double **rdfs, double **rdfs_d, tW_RDF_selection *sel, int n_rdfs, double bw, int n_bins, bool derv)
{
  int ii, i, jj, j, k, itype, jtype, idx, rdf_idx;
  double dist;
  double volume = fr->contents->box[0][0] * fr->contents->box[1][1] * fr->contents->box[2][2];
  double factor, dV, r;
  rdf_idx = 0;
  for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
  {
    itype = sel[rdf_idx].typeA;
    jtype = sel[rdf_idx].typeB;
    if (itype == jtype)
    {
      int hold = top->idx_file->n_at_per_group[itype];
      factor = volume / (double)((double)(hold * hold) / 2.0);
      if (fr->counter == 1) { fprintf(stderr,"factor_ii: %g \n",factor); }
    }
    else
    {
      factor = volume / (double)(top->idx_file->n_at_per_group[itype] * top->idx_file->n_at_per_group[jtype]);
      if (fr->counter == 1) { fprintf(stderr,"factor_ij: %g \n",factor); }
    }
    rdfs[rdf_idx][0] = 0.0;
    for ( k = 1; k < n_bins; ++k)
    {
      r = (double) k * bw;
      dV = ((double) (4.0 * M_PI * r * r * bw));
      rdfs[rdf_idx][k] *= factor / dV;
      if (derv) { rdfs_d[rdf_idx][k] *= factor / dV; }
    }
  }
}

void do_rdfs(tW_gmx_topology *top, tW_gmx_trxframe *fr, double **rdfs, double **rdfs_d, tW_RDF_selection *sel, int n_rdfs, double bw, int n_bins, bool scale_now, bool derv, double ener)
{
  int ii, i, jj, j, k, itype, jtype, idx, rdf_idx;
  double dist;
  dvec dr;
  for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
  {
    itype = sel[rdf_idx].typeA;
    jtype = sel[rdf_idx].typeB;
    for (ii = 0; ii < top->idx_file->n_at_per_group[itype]; ++ii)
    {
      i = top->idx_file->atoms_in_group[itype][ii];
      printf("This is the structure stored in top->idx_file_atoms_in_group[itype][ii] \n");
      printf("n_at_per_group[itype]: %d atoms_in_group[itype][ii]: %d \n", top->idx_file->n_at_per_group[itype], top->idx_file->atoms_in_group[itype][ii]); 
      if (itype == jtype)
      {
        for (jj = ii+1; jj < top->idx_file->n_at_per_group[jtype]; ++jj)
        {
          j = top->idx_file->atoms_in_group[jtype][jj];
          dist = get_rect_prism_dx(fr->contents->x[i], fr->contents->x[j], fr->contents->box, dr);
          idx = get_idx(dist, bw, n_bins);
          if (idx != -1) 
          { 
            rdfs[rdf_idx][idx] += 1.0;
            if (derv)
            {
              rdfs_d[rdf_idx][idx] += ener;
            }
          }
        }
      }
      else
      {
        for (jj = 0; jj < top->idx_file->n_at_per_group[jtype]; ++jj)
        {
          j = top->idx_file->atoms_in_group[jtype][jj];
          dist = get_rect_prism_dx(fr->contents->x[i], fr->contents->x[j], fr->contents->box, dr);
          idx = get_idx(dist, bw, n_bins);
          if (idx != -1)
          { 
            rdfs[rdf_idx][idx] += 1.0;
            if (derv)
            {
            rdfs_d[rdf_idx][idx] += ener;
            }
          }
        }
      }
    }
  }

/* 
I think we have to do the scaling on a per frame basis incase it's an NPT simulation...
If it was NVT we could save it all for the end, which would be much more computationally efficient.
But if it's not...

Let's make it a command line arg?
 
*/
  if (scale_now)
  {
    scale_rdfs(top, fr, rdfs, rdfs_d, sel, n_rdfs, bw, n_bins, derv);
    /*double volume = fr->contents->box[0][0] * fr->contents->box[1][1] * fr->contents->box[2][2];
    double factor, dV, r;
    rdf_idx = 0;
    for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
    {
      itype = sel[rdf_idx].typeA;
      jtype = sel[rdf_idx].typeB;
      if (itype == jtype) 
      { 
        int hold = top->idx_file->n_at_per_group[itype];
        factor = volume / (double)(hold * (hold) / 2); 
        if (fr->counter == 1) { fprintf(stderr,"factor_ii: %g \n",factor); }
      }
      else
      {
        factor = volume / (double)(top->idx_file->n_at_per_group[itype] * top->idx_file->n_at_per_group[jtype]);
        if (fr->counter == 1) { fprintf(stderr,"factor_ij: %g \n",factor); }
      }
      rdfs[rdf_idx][0] = 0.0;
      for ( k = 1; k < n_bins; ++k)
      {
        r = (double) k * bw;
        dV = ((double) (4.0 * M_PI * r * r * bw));
        rdfs[rdf_idx][k] *= factor / dV;
      }
    }*/
  }
}

void increment_total_rdfs(double **total_rdfs, double **rdfs, double **total_rdfs_d, double **rdfs_d, int n_frames_local, int n_rdfs, int n_bins, bool derv)
{
  int i, j;
  double factor = ((double) (n_frames_local - 1) / (double) (n_frames_local));
  for (i = 0; i < n_rdfs; ++i)
  {
    for (j = 0; j < n_bins; ++j)
    {
      total_rdfs[i][j] = total_rdfs[i][j] * factor + rdfs[i][j] / n_frames_local;
      rdfs[i][j] = 0.0;

      if (derv)
      {
      total_rdfs_d[i][j] = total_rdfs_d[i][j] * factor + rdfs_d[i][j] / n_frames_local;
      rdfs_d[i][j] = 0.0;
      }
    }
  }
}

void save_rdfs(double **rdfs, tW_gmx_topology *top, tW_RDF_selection *sel, int n_rdfs, int n_bins, double bw, int frame_nr, bool ttl)
{
  int i, j, k, rdf_idx;
  double r;
  for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
  {
    i = sel[rdf_idx].typeA;
    j = sel[rdf_idx].typeB;
    char *fnm = (char *) ecalloc(100,sizeof(char));
    if (ttl) { sprintf(fnm,"BocsRDF_%s_%s_total.dat",top->idx_file->group_names[i],top->idx_file->group_names[j]); }
    else { sprintf(fnm,"BocsRDF_%s_%s_fr_%d.dat",top->idx_file->group_names[i],top->idx_file->group_names[j],frame_nr); }
    FILE *fp = open_file(fnm,'w');
    for (k = 0; k < n_bins; ++k)
    {
      r = (double) k * bw;
      fprintf(fp,"%f %g \n",r,rdfs[rdf_idx][k]);
    }     
    fclose(fp);
    efree(fnm);
  }
}

void save_rdfs_d(double **rdfs_d, tW_gmx_topology *top, tW_RDF_selection *sel, int n_rdfs, int n_bins, double bw, int frame_nr, bool ttl)
{
  int i, j, k, rdf_idx;
  double r;
  for (rdf_idx = 0; rdf_idx < n_rdfs; ++rdf_idx)
  {
    i = sel[rdf_idx].typeA;
    j = sel[rdf_idx].typeB;
    char *fnm = (char *) ecalloc(100,sizeof(char));
    if (ttl) { sprintf(fnm,"BocsRDF_D_%s_%s_total.dat",top->idx_file->group_names[i],top->idx_file->group_names[j]); }
    else { sprintf(fnm,"BocsRDF_D_%s_%s_fr_%d.dat",top->idx_file->group_names[i],top->idx_file->group_names[j],frame_nr); }
    FILE *fp = open_file(fnm,'w');
    for (k = 0; k < n_bins; ++k)
    {
      r = (double) k * bw;
      fprintf(fp,"%f %g \n",r,rdfs_d[rdf_idx][k]);
    }
    fclose(fp);
    efree(fnm);
  }
}
