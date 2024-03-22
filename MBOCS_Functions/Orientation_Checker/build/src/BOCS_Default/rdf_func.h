#ifndef RDF_FUNC_H
#define RDF_FUNC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"

#include "safe_mem.h" 
#include "io_read.h"
#include "gromacs_topology.h"
#include "LDD.h"

typedef struct tW_RDF_selection
{
  int typeA;
  int typeB;
} tW_RDF_selection;

void get_rdf_selections(tW_gmx_topology *top,tW_RDF_selection *sel, int n_rdfs);

void scale_rdfs(tW_gmx_topology *top, tW_gmx_trxframe *fr, double **rdfs, double **rdfs_d, tW_RDF_selection *sel, int n_rdfs, double bw, int n_bins, bool derv);

void do_rdfs(tW_gmx_topology *top, tW_gmx_trxframe *fr, double **rdfs, double **rdfs_d, tW_RDF_selection *sel, int n_rdfs, double bw, int n_bins, bool scale_now, bool derv, double ener);

void increment_total_rdfs(double **total_rdfs, double **rdfs, double **total_rdfs_d, double **rdfs_d, int n_frames_local, int n_rdfs, int n_bins, bool derv);

void save_rdfs(double **rdfs, tW_gmx_topology *top, tW_RDF_selection *sel, int n_rdfs, int n_bins, double bw, int frame_nr, bool ttl);

void save_rdfs_d(double **rdfs_d, tW_gmx_topology *top, tW_RDF_selection *sel, int n_rdfs, int n_bins, double bw, int frame_fr, bool ttl);

#ifdef __cplusplus
}
#endif

#endif
