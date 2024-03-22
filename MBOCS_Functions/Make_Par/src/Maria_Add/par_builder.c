// par_builder is a functionality constructed by MCL to make forcematching larger molecules easier
// In this early draft, this script works off of the BOCS.btp file so that we can incorporate non-gromacs MD packages into it using the translator later. 
// I haven't included LD capability in this script and haven't implemented Michael's patch where molblocks != moltypes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../BOCS_Default/cgff_types.h"
#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/gromacs_topology.h"
#include "par_builder_fn.h"
#include "../BOCS_Default/wnoid_math.h"

int main(int argc, char * argv[])
{

  /*SET UP USER INPUT*/
  tW_word fnmi, fnmo; // Strings to store the names of files [in and out]

  tW_gmx_trxframe *fri = init_tW_gmx_trxframe(); // These are data structures we've made room for in memory, we may or may not fill them. 
  tW_gmx_trxframe *fro = init_tW_gmx_trxframe();
  tW_gmx_topology *top = init_tW_gmx_topology();

  int n_params = 2; // Expect the user to pass to 2 args when using this function
  int n_arg_found; // Stores how many args we actually get
  int i = 0;
  t_pargs *params = (t_pargs *) ecalloc(n_params, sizeof(t_pargs)); // Make room for the three args we expect together with handy members for tracking reqs/ if they've been found

  init_arg(&(params[0]),"-f",etSTRING,"input filename [.btp file] (req.)",TRUE); // This is where we specify how those handy structures get filled
  init_arg(&(params[1]),"-o",etSTRING,"output filename (req.)",TRUE);

  n_arg_found = get_command_line_args(argc, argv, n_params, params); // This extracts how many arguments we got from the user and fills in their values, we set this to -9999 when -h is passed

  print_arg_table(n_params,params,(n_arg_found > 0 ? FALSE : TRUE)); // This prints what arguments we have and their current values / the associated warning notes  

  if (n_arg_found <= 0) { return 0; } // This breaks our program if the user asked for help with -h

  check_mand_args(params,n_params);

  /*Now We'll get specific to this function*/
//  if (strstr(params[1].value,"par.txt") != NULL) // check if they get that this makes par.txt
//  {
//    fprintf(stderr,"Warning: cgff will require par.txt as input, this script assumes you want par.txt from a .btp file.\n");

    top->eFileType = eBOCS; // making sure read_topology gets the right type 
    top->b_tpr = read_topology(top,params[0].value); // Wrapper function, fills top from input based on eFileType, here it calls read_tpr_dump, returns true if reading was successful, Michael wrote this
  
    if (top->b_tpr) { write_par(params[1].value,top); }
    else
    {
      fprintf(stderr,"ERROR: there was a problem trying to read the dumped topology file %s \n",params[0].value);
      return 1;
    }
 //   return 0;
//  }


  return 0;
}
