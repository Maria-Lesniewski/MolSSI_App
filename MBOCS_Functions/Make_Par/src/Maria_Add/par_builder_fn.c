#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../BOCS_Default/cgff_types.h"
#include "../BOCS_Default/gromacs_topology.h"
#include "../BOCS_Default/safe_mem.h"
#include "../BOCS_Default/io_read.h"
#include "../BOCS_Default/wnoid_math.h"
#include "par_builder_fn.h"
// It helps to read this file bottom -> top [top are the last functions I wrote]

/******************************************************************************************
 * print_imnb(file, pair_inter_list, list_len) : prints inter_name and involved atoms for 
 *                                               Intramolecular nb_pair interactions 
 ******************************************************************************************/
void print_imnb( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "[IntraMolec_NB_Pair] %d\n", len);
  fprintf(fp, "%16s %8s %8s\n", "! inter_name",   "type1",  "type2");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp, "Intra_%s%s %8s %8s\n", ilist[i][0], ilist[i][1], ilist[i][0], ilist[i][1]);
  }
  fprintf(fp, "[End IntraMolec_NB_Pair]\n");
}

/*****************************************************************************************
 * print_dihedral(file, pair_inter_list, list_len) : prints inter_name and involved atoms for 
 *                                                   dihedral interactions in  the par file
 ******************************************************************************************/
void print_dihedral( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "[Dihedral] %d\n", len);
  fprintf(fp, "%16s %8s %8s %8s\n", "! inter_name",   "type1",  "type2", "type3");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp, "Dihedral_%s-%s-%s-%s %8s %8s %8s %8s\n", ilist[i][0], ilist[i][1], ilist[i][2], ilist[i][3], ilist[i][0], ilist[i][1], ilist[i][2], ilist[i][3]);
  }
  fprintf(fp, "[End Dihedral]\n");
}

/*****************************************************************************************
 * print_angle(file, pair_inter_list, list_len) : prints inter_name and involved atoms for 
 *                                                angular interactions in  the par file
 ******************************************************************************************/
void print_angle( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "[Angle] %d\n", len);
  fprintf(fp, "%16s %8s %8s %8s\n", "! inter_name",   "type1",  "type2", "type3");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp, "Angle_%s-%s-%s %8s %8s %8s\n", ilist[i][0], ilist[i][1], ilist[i][2], ilist[i][0], ilist[i][1], ilist[i][2]);
  }
  fprintf(fp, "[End Angle]\n");
}

/******************************************************************************************
 * print_bond(file, pair_inter_list, list_len) : prints inter_name and involved atoms for 
 *                                               bonded interactions in  the par file
 ******************************************************************************************/
void print_bond( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "[BondStretch] %d\n", len);
  fprintf(fp, "%16s %8s %8s\n", "! inter_name",   "type1",  "type2");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp, "Bond_%s-%s %8s %8s\n", ilist[i][0], ilist[i][1], ilist[i][0], ilist[i][1]);
  }
  fprintf(fp, "[End BondStretch]\n");
}

/******************************************************************************************
 * print_pair(file, pair_inter_list, list_len) : prints inter_name and involved atoms for 
 *                                               nb_pair interactions in  the par file
 ******************************************************************************************/
void print_pair( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "[Pair_Interaction] %d\n", len);
  fprintf(fp, "%16s %8s %8s\n", "! inter_name",   "type1",  "type2");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp, "%8s%-8s %8s %8s\n", ilist[i][0], ilist[i][1], ilist[i][0], ilist[i][1]);
  }
  fprintf(fp, "[End Pair_Interaction]\n");
}

void write_interactions(FILE *fp, tW_word** ilist, int len, int mode)
{
  switch(mode)
  {
    case mNB:
    {
      print_pair(fp, ilist,len);
      break;
    }
    case mBOND:
    {
      print_bond(fp, ilist, len);
      break;
    }
    case mANGLE:
    {
      print_angle(fp, ilist, len);
      break;
    }
    case mDIH:
    {
      print_dihedral(fp, ilist, len);
      break;
    }
    case mIMNB:
    {
      print_imnb(fp, ilist,len);
      break;
    }
  }
}

/**********************************************************************************************************************
 * total_print_imnb(file, impropernb_list, num_interactions) : prints the list of improper nb interactions and the 
 *                                                             default Bspline basis settings for forcematching 
 *********************************************************************************************************************/
void total_print_imnb( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "%16s %18s %8s %8s %8s %8s %8s %12s \n", "!inter_name","inter_type", "basis", "dr", "R_min", "R_max", "spline_order", "n_smooth");
  printf("SELECTING DEFAULT IMNB BEHAVIOR CUBIC BSPLINE r = 0.00 - 2.00 bw=0.02 \n");
  printf("Total length: %d\n", len);
  for (int i = 0; i < len; i++)
  {
    fprintf(fp,"Intra_%s%s %18s %8s %8s %8s %8s %8s %8s\n", ilist[i][0], ilist[i][1], "IntraMolec_NB_Pair", "Bspline", "0.02", "0.000", "2.000", "4", "0");
  }
}
/**********************************************************************************************************************
 * total_print_dihedral(file, dihedral_list, num_interactions) : prints the list of dihedral interactions and the 
 *                                                                 default linear basis settings for forcematching 
 *********************************************************************************************************************/
void total_print_dihedral( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "!inter_name %10s %8s %8s %8s %8s %12s \n", "inter_type", "basis", "dr", "R_min", "R_max", "n_smooth");
  printf("SELECTING DEFAULT DIHEDRAL BEHAVIOR linear r = -180 to 180 bw=0.500 \n");
  for (int i = 0; i < len; i++)
  {
    fprintf(fp,"Dihedral_%s-%s-%s-%s %14s %8s %8s %8s %8s %8s \n", ilist[i][0], ilist[i][1], ilist[i][2], ilist[i][3], "Dihedral", "linear", "0.500", "-180", "180", "0");
  }
}
/**********************************************************************************************************************
 * total_print_angle(file, angular_list, num_interactions) : prints the list of angular interactions and the 
 *                                                           default linear basis settings for forcematching 
 *********************************************************************************************************************/
void total_print_angle( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "!inter_name %10s %8s %8s %8s %8s %12s \n", "inter_type", "basis", "dr", "R_min", "R_max", "n_smooth");
  printf("SELECTING DEFAULT ANGLE BEHAVIOR linear r = 0.00 - 180 bw=0.500 \n");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp,"Angle_%s-%s-%s %14s %8s %8s %8s %8s %8s \n", ilist[i][0], ilist[i][1], ilist[i][2], "Angle", "linear", "0.500", "0.000", "180", "0");
  }
}
/**********************************************************************************************************************
 * total_print_bond(file, angular_list, num_interactions) : prints the list of bonded interactions and the 
 *                                                           default linear basis settings for forcematching 
 *********************************************************************************************************************/
void total_print_bond( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "!inter_name %10s %8s %8s %8s %8s %12s \n", "inter_type", "basis", "dr", "R_min", "R_max", "n_smooth");
  printf("SELECTING DEFAULT BOND BEHAVIOR linear r = 0.000 - 0.7 bw=0.002 \n");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp,"Bond_%s-%s %14s %8s %8s %8s %8s %8s \n", ilist[i][0], ilist[i][1], "BondStretch", "linear", "0.002", "0.000", "0.7",  "0");
  }
}
/**********************************************************************************************************************
 * total_print_pair(file, angular_list, num_interactions) : prints the list of nonbonded pair interactions and the 
 *                                                           default Bspline basis settings for forcematching 
 *********************************************************************************************************************/
void total_print_pair( FILE *fp, tW_word** ilist, int len)
{
  fprintf(fp, "%16s %16s %8s %8s %8s %8s %8s %12s \n", "!inter_name","inter_type", "basis", "dr", "R_min", "R_max", "spline_order", "n_smooth");
  printf("SELECTING DEFAULT NB BEHAVIOR CUBIC BSPLINE r = 0.00 - 2.00 bw=0.02 \n");
  for (int i = 0; i < len; i++)
  {
     fprintf(fp,"%8s%-8s %16s %8s %8s %8s %8s %8s %8s\n", ilist[i][0], ilist[i][1], "Pair_Interaction", "Bspline", "0.02", "0.000", "2.000", "4", "0");
  }
}
/**********************************************************************************************************************
 * write_total_list(file, angular_list, num_interactions) : wrapper function for printing all interaction type options and 
 *                                                          default basis settings. [change this function as more types added]
 *********************************************************************************************************************/
void write_total_list( FILE *fp, tW_word** ilist, int len, int mode)
{
  switch(mode)
  {
    case mNB:
    {
      total_print_pair(fp, ilist, len); 
      break;
    }
    case mBOND:
    {
      total_print_bond(fp, ilist, len);
      break;
    }
    case mANGLE:
    {
      total_print_angle(fp, ilist, len);
      break;
    }
    case mDIH:
    {
      total_print_dihedral(fp, ilist, len);
      break;
    }
    case mIMNB:
    {
      total_print_imnb(fp, ilist, len);
      break;
    }
    default:
    {
      fprintf(fp, "[End Inter_Types]");
    }
  }
}

/********************************************************************************
 * fill_impropernb(inter_list, my_moltype, found_vector, index_to_check) : 
 * Analogous to search_improper, but for list filling instead of counting
 * loops through all the angles between atoms in a molecule, determines if the
 * types are unique and if so adds the three types to the angular inter_list
 *********************************************************************************/
void fill_impropernb(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check)
{
  //Find the type of atoms involved in this pair
  tW_word typeA, typeB;

  strcpy(typeA,  mol->atom_types[mol->imnb_ij[check][0]]);
  strcpy(typeB, mol->atom_types[mol->imnb_ij[check][1]]);

  //See if any other improper listed involve these same two types
  for (int i = 0; i < mol->n_imnbs; i++)
  {
    if (
        strcmp(typeA, mol->atom_types[mol->imnb_ij[i][0]]) == 0 && strcmp(typeB, mol->atom_types[mol->imnb_ij[i][1]]) == 0 ||
        strcmp(typeA, mol->atom_types[mol->imnb_ij[i][1]]) == 0 && strcmp(typeB, mol->atom_types[mol->imnb_ij[i][0]]) == 0
       )
       {
         beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
       }
  }
  strcpy(ilist[current_index][0], typeA);// Take the site types from this instance
  strcpy(ilist[current_index][1], typeB);  
}

void fill_impropernb_list(tW_gmx_topology* top, tW_word **ilist, int list_len)
{
  int current_index = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the pair before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_imnbs, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_imnbs; k++){ repeat_found[k] = FALSE;}

    //Search through pairs, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_imnbs; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        fill_impropernb(ilist, current_index, &(top->molecules[i]), repeat_found, j);
        current_index++;
      }
    }
  }
  printf("Finished filling improper nb interactions, found %d interactions\n", current_index);
}

/********************************************************************************
 * fill_angles(inter_list, my_moltype, found_vector, index_to_check) : 
 * Analogous to search_dihedrals, but for list filling instead of counting
 * loops through all the dihedrals between atoms in a molecule, determines if the
 * types are unique and if so adds the four types to the dih inter_list
 *********************************************************************************/
void fill_dihedrals(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check)
{
  //Find the type of atoms involved in this dihedral
  tW_word nameA, nameB, nameC, nameD;
  strcpy(nameA,  mol->atom_types[mol->dih_ijkl[check][0]]);
  strcpy(nameB, mol->atom_types[mol->dih_ijkl[check][1]]); // Central
  strcpy(nameC, mol->atom_types[mol->dih_ijkl[check][2]]); // Pair Types
  strcpy(nameD, mol->atom_types[mol->dih_ijkl[check][3]]);
  //See if any other dihedral instances involve these same types
  for (int i = 0; i < mol->n_dihs; i++)
  {
   /* if (
        (strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 && strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 )||
        (strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 && strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 )
       )
       {
         if (
            (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][0]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][3]]) == 0) ||
            (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][3]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][0]]) == 0)
            ) need to put back a }*/
     if ( 
         (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][0]]) == 0 && strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 &&
          strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][3]]) == 0) ||
         (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][3]]) == 0 && strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 &&
          strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][0]]) == 0)
        )
         {
              beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
         }
  }
  strcpy(ilist[current_index][0], nameA);// Take the site types from this instance
  strcpy(ilist[current_index][1], nameB);
  strcpy(ilist[current_index][2], nameC);
  strcpy(ilist[current_index][3], nameD); // need to update analogous functions so that not contanstly re-initializing same item 
}
/***************************************************************************************
 * fill_dihedral_list(system_top, inter_list, length_inter_list) : Analogous to get_n_dihedrals
 *                                                              but fills unique interactions
 *                                                              found into inter_list passed
 ***************************************************************************************/
void fill_dihedral_list(tW_gmx_topology* top, tW_word **ilist, int list_len)
{
  int current_index = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the dihedral before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_dihs, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_dihs; k++){ repeat_found[k] = FALSE;}

    //Search through dihedrals, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_dihs; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        fill_dihedrals(ilist, current_index, &(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
        current_index++;
      }
    }
  }
  printf("Finished searching for dihedrals, found %d interactions\n", current_index);
}

/********************************************************************************
 * fill_angles(inter_list, my_moltype, found_vector, index_to_check) : 
 * Analogous to search_angles, but for list filling instead of counting
 * loops through all the angles between atoms in a molecule, determines if the
 * types are unique and if so adds the three types to the angular inter_list
 *********************************************************************************/
void fill_angles(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check)
{
  //Find the type of atoms involved in this angle
  tW_word nameA, nameB, nameC;
  strcpy(nameA,  mol->atom_types[mol->angle_ijk[check][0]]);
  strcpy(nameB, mol->atom_types[mol->angle_ijk[check][1]]); // Central Atom
  strcpy(nameC, mol->atom_types[mol->angle_ijk[check][2]]);
  //See if any other instance of same angle type is present
  for (int i = 0; i < mol->n_angles; i++)
  {
    if (strcmp(nameB, mol->atom_types[mol->angle_ijk[i][1]]) == 0)
    {
      if (
          (strcmp(nameA, mol->atom_types[mol->angle_ijk[i][0]]) == 0 && strcmp(nameC, mol->atom_types[mol->angle_ijk[i][2]]) == 0) ||
          (strcmp(nameA, mol->atom_types[mol->angle_ijk[i][2]]) == 0 && strcmp(nameC, mol->atom_types[mol->angle_ijk[i][0]]) == 0)
         )
         {
           beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
         }
    }
  }
  strcpy(ilist[current_index][0], nameA);// Take the site types from this instance
  strcpy(ilist[current_index][1], nameB);
  strcpy(ilist[current_index][2], nameC); 
}

/***************************************************************************************
 * fill_angle_list(system_top, inter_list, length_inter_list) : Analogous to get_n_angles
 *                                                              but fills unique interactions
 *                                                              found into inter_list passed
 ***************************************************************************************/
void fill_angle_list(tW_gmx_topology* top, tW_word **ilist, int list_len)
{
  int current_index = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the angle before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_angles, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_angles; k++){ repeat_found[k] = FALSE;}

    //Search through angles, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_angles; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        fill_angles(ilist, current_index, &(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
        current_index++;
      }
    }
  }
  printf("Finished searching for angles, found %d interactions\n", current_index);
}


/**********************************************************************************************
 * fill_bonds(bond_inter_list, bonds_found_so_far, my_molecule_type, duplicate_tracker, bond_instance_index_to_check) :
 * Analogous to search_bonds, but fills in the site_types found in the passed list instead of just counting unique pairs
 **************************************************************************************************/
void fill_bonds(tW_word **ilist, int current_index, tW_molecule* mol, bool *beq_found, int check)
{
  //Find names we need to check
  char nameA[50], nameB[50];
  strcpy(nameA,  mol->atom_types[mol->bond_ij[check][0]]);
  strcpy(nameB, mol->atom_types[mol->bond_ij[check][1]]);

  //See if any other bonds involve these same two types
  for (int i = 0; i < mol->n_bonds; ++i)
  {
    if (
       (strcmp(nameA, mol->atom_types[mol->bond_ij[i][0]]) == 0 && strcmp(nameB, mol->atom_types[mol->bond_ij[i][1]]) == 0 )
       ||
       (strcmp(nameA, mol->atom_types[mol->bond_ij[i][1]]) == 0 && strcmp(nameB, mol->atom_types[mol->bond_ij[i][0]]) == 0)
       )
       {
	 beq_found[i]=TRUE; // Mark all duplicates
       } 
   }
  strcpy(ilist[current_index][0], nameA);
  strcpy(ilist[current_index][1], nameB); // Take types from this instance
}

/**********************************************************************************************
 * fill_bond_list(system_top, bond_inter_list, goal_length) : Analogous to get_n_distinct_bonds
 *                                                            but fills the interaction types into 
 *                                                            inter_list instead of counting
 *********************************************************************************************/
void fill_bond_list(tW_gmx_topology* top, tW_word **ilist, int list_len)
{
  int current_index = 0;

  for (int i = 0; i < top->n_molecule_types; i++) // check bond total for each molecule type
  {
    // Store if we've checked a pair in the molecule before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_bonds, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_bonds; k++){ repeat_found[k] = FALSE;}

    // Search pairs, not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_bonds; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        fill_bonds(ilist, current_index, &(top->molecules[i]), repeat_found, j); // sets repeat_found elements to true
	current_index++;
      }
    }
  }
  printf("Finished bond search, found %d bonded interactions\n", current_index);
}
/**********************************************************************************
 *fill_intermol_pairs(pair_inter_list, my_moltype, num_inter_found_already) :
 * Finds pair interactions from different types of molecules interacting,
 * Fills inter_list and returns the number of inter found [filled]
 **********************************************************************************/
int fill_intermol_pairs(tW_word **ilist, tW_gmx_topology* top, int ifoundsofar)
{
  int index = ifoundsofar;
  int count_added = 0;

  for (int i = 0; i < top->n_molecule_types; i++)
  {
    for (int j = i + 1; j < top->n_molecule_types; j++) // Scan through molecule types
    {
      for (int k = 0; k < top->molecules[i].n_apm; k++) // Scan each atom on each molecule for types
      {
        for (int l = 0; l < top->molecules[j].n_apm; l++)
        {
          if (check_pairlist(top->molecules[i].atom_types[k], top->molecules[j].atom_types[l], ilist, ifoundsofar) == FALSE)
          {
            strcpy(ilist[index][0], top->molecules[i].atom_types[k]);
            strcpy(ilist[index][1], top->molecules[j].atom_types[l]);
            
            index++;
            count_added++;
          }
        }
      }
    }
  }
  return count_added;
}

/*************************************************************************************
 *molecule_choose_2_interactions(pair_inter_list, my_moltype, num_inter_found_already) :
 * Finds pair interactions associated with molecules that see another copy of themselves
 * inside of the simulation [naive listing works], fills list & returns num inter found
 *************************************************************************************/
int molecule_choose_2_interactions(tW_word **ilist, tW_molecule *mol, int ifoundsofar)
{
  int index = ifoundsofar;
  int count_added = 0;
  
  for (int i = 0; i < mol->n_apm; i++)
  {
    for (int j = i; j < mol->n_apm; j++)
    {
      if (check_pairlist(mol->atom_types[i], mol->atom_types[j], ilist, ifoundsofar) == FALSE)
      {
        strcpy(ilist[index][0], mol->atom_types[i]);
        strcpy(ilist[index][1], mol->atom_types[j]);
        printf("%s %s  ADDED \n", mol->atom_types[i], mol->atom_types[j]);
	printf("ADDING TO COUNT %d -> %d\n", count_added, count_added+1);
        index++;
        count_added++;
      }
    }
  }
  printf("Count returned  %d\n", count_added); 
  return count_added;  
}

/*****************************************************************************************
 *fill_nb_list_unique_mol(interaction_list, my_molecule_type, num_inter_found_already) :
 * Finds pair interactions associated with long molecules that there are only
 * one of in the simulation [accounts for exclusions], fills list & returns num inter found
 ****************************************************************************************/
int fill_nb_list_unique_mol(tW_word **ilist, tW_molecule *mol, int ifoundsofar)
{
  int index = ifoundsofar;
  int count_added = 0;

  // look at all the atoms/sites in the molecule
  for (int i = 0; i < mol->n_apm; i++) // sit on an atom
  {
    for (int j = mol->n_apm-1; j>i; j--) // look at all the other atoms in this molecule we havent sat on yet 
    {
      if (check_list(mol->excls[i],mol->n_epa[i],j) == FALSE) // If we need to consider this atom pair
      {
        // check if we have the types in ilist alread
        if (check_pairlist(mol->atom_types[i], mol->atom_types[j], ilist, ifoundsofar) == FALSE)
        {
          // We need to add these atom_types to the list and make room for more
          strcpy(ilist[ifoundsofar][0], mol->atom_types[i]);
          strcpy(ilist[ifoundsofar][1], mol->atom_types[j]);
      
      	  ifoundsofar++;
          count_added++;
        }
      }
    }
  }
  return count_added; 
}

/**********************************************************************************************
 * brutal_list_pair_int(system_top, interaction_list, goal_length) : lists the required pair int
 *                                                                   in the case that exclusions
 *                                                                   need to be considered 
 **********************************************************************************************/
void brutal_list_pair_int(tW_gmx_topology *top, tW_word **ilist, int list_len)
{
  //Seperate how we count interactions in molecules types with nmol > 1 from nmol = 1 types
  int* lone_indices = (int*) emalloc(sizeof(int));
  int* multi_indices = (int*) emalloc(sizeof(int));
  int list_lengths[2];
  sort_moltypes(top, lone_indices, multi_indices, list_lengths);
  
  //List the interactions for site types within a molecule with nmol = 1
  int interactions_found = 0;
  for (int i = 0; i < (list_lengths[0]-1); i++)
  {
   interactions_found = interactions_found + fill_nb_list_unique_mol(ilist, &(top->molecules[lone_indices[i]]), interactions_found);
  }
 
  //List the interactions for site types contained in molecules with nmol > 1
  for (int j = 0; j < (list_lengths[1]-1); j++)
  {
   interactions_found = interactions_found + molecule_choose_2_interactions(ilist, &(top->molecules[multi_indices[j]]), interactions_found);
  }

  //List the cross interactions for site types between different molecules
  interactions_found = interactions_found + fill_intermol_pairs(ilist, top, interactions_found);

  printf("Finished nb search, found %d pairs\n", interactions_found);
}

void list_choose_2_interactions(tW_gmx_topology *top, tW_word **ilist, int list_len)
{
  int current_index = 0;
  // A-B interactions on j!=i, A-A on j = i
  for (int i = 0; i < top->n_atomtypes; i++)
  {
    for (int j = i; j < top->n_atomtypes; j++)
    {
      strcpy(ilist[current_index][0], top->atom_type_names[i]);
      strcpy(ilist[current_index][1], top->atom_type_names[j]);
      current_index++;
    } 
  }
  //printf("END OF LISTING FILLED %d pairs\n", current_index);
}

/***********************************************************************************************************
 * fill_nb_list(system_top, pair_interaction_list, goal_length) : Searches for pair nonbonded interactions 
 *                                                                for an initialized list of known length
 *                                                                Fills pair names into list
 ***********************************************************************************************************/
void fill_nb_list(tW_gmx_topology *top, tW_word **ilist, int list_len)
{
  // Different algorithms apply for listing interactions 
  if (naive_exclusion_check(top) == TRUE)
  {
     list_choose_2_interactions(top, ilist, list_len);
  }
  else
  {
     brutal_list_pair_int(top, ilist, list_len);
  } 
}
/******************************************************************************************************************
 * search_improper( my_mol, duplicate_found_vector, interaction_instance_i'm_checking) : loops through all instances
 *                                                                                       of imnb interactions,
 *                                                                                       sets beq to true where eq.
 *                                                                                       types have been found
 ******************************************************************************************************************/
int search_improper(tW_molecule * mol, bool * beq_found, int check)
{
  //Find the type of atoms involved in this pair
  tW_word typeA, typeB;
  
  strcpy(typeA,  mol->atom_types[mol->imnb_ij[check][0]]);
  strcpy(typeB, mol->atom_types[mol->imnb_ij[check][1]]);
  
  //See if any other improper listed involve these same two types
  for (int i = 0; i < mol->n_imnbs; i++)
  {
    if (
        strcmp(typeA, mol->atom_types[mol->imnb_ij[i][0]]) == 0 && strcmp(typeB, mol->atom_types[mol->imnb_ij[i][1]]) == 0 ||
        strcmp(typeA, mol->atom_types[mol->imnb_ij[i][1]]) == 0 && strcmp(typeB, mol->atom_types[mol->imnb_ij[i][0]]) == 0
       )
       {
         beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
       }
  }
  return 1;
}

/************************************************************************************************
 * get_n_impropernb(system_top) : returns the number of improper nonbonded (e.g. LJ14) interactions
 *                                required by the input topology
 ************************************************************************************************/
int get_n_impropernb(tW_gmx_topology *top)
{
  int total_improper = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the pair before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_imnbs, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_imnbs; k++){ repeat_found[k] = FALSE;}

    //Search through pairs, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_imnbs; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        total_improper = total_improper + search_improper(&(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
      }
    }
  }
  return total_improper;
}


/**************************************************************************************************************
 * list_intermol_interactions(interaction_list, system_top, length_interact_list) : Adds A-B interactions to 
 *										    list of interactions,
 *										    returns # interactions
 *										    appended.
 *										    Here A is an atomtype in moltype1
 *										         B is an atomtype in moltype2
 **************************************************************************************************************/
int list_intermol_interactions(tW_word ** ilist, tW_gmx_topology *top, int* len_ilist)
{
  int length_added = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    printf("Checking inter neighbors for moltype: %s\n", top->molecules[i].molname);
    for (int j = i + 1; j < top->n_molecule_types; j++) // Scan through molecule types
    {
      printf("Its molecule buddy is %s \n", top->molecules[j].molname);
      for (int k = 0; k < top->molecules[i].n_apm; k++) // Scan each atom on each molecule for types
      {
        for (int l = 0; l < top->molecules[j].n_apm; l++)
	{
	  if (check_pairlist(top->molecules[i].atom_types[k], top->molecules[j].atom_types[l], ilist, *len_ilist) == FALSE)
	  {
	    strcpy(ilist[*len_ilist-1][0], top->molecules[i].atom_types[k]);
            strcpy(ilist[*len_ilist-1][1], top->molecules[j].atom_types[l]);

            *len_ilist= *len_ilist + 1;
            printf("Allocating Room for Entry %d\n", *len_ilist);
            ilist = (tW_word **) erealloc(ilist, *len_ilist * sizeof(tW_word));
            ilist[*len_ilist-1] = (tW_word*) erealloc(ilist[*len_ilist-1], 2 * sizeof(tW_word));
            length_added++;
	  }
	}
      }
    }
  }
  return length_added;
}


/*****************************************************************************************************
 * list_cross_interactions(listofint, my_mol, len_list) : Adds A-B interactions between molecules of 
 *                                                        same type to the list of pair interactions.
 *                                                        [A & B are atomtypes within the molecule]
 *                                                        returns the number of interactions appended
 *                                                        to the list
 *****************************************************************************************************/
int list_cross_interactions(tW_word ** ilist, tW_molecule * mol, int* len_ilist)
{
  int length_added = 0;
  for (int i = 0; i < mol->n_apm; i++)
  {
    for (int j = i+1; j < mol->n_apm; j++)
    {
      if (check_pairlist(mol->atom_types[i], mol->atom_types[j], ilist, *len_ilist) == FALSE)
      {
        strcpy(ilist[*len_ilist-1][0], mol->atom_types[i]);
        strcpy(ilist[*len_ilist-1][1], mol->atom_types[j]);

        *len_ilist= *len_ilist + 1;
        ilist = (tW_word **) erealloc(ilist, *len_ilist * sizeof(tW_word));
        ilist[*len_ilist-1] = (tW_word*) erealloc(ilist[*len_ilist-1], 2 * sizeof(tW_word));
        length_added++;
      }
    }
  }
  return length_added;
}

/*****************************************************************************************************
 * list_self_interactions(list_of_interactions, my_molecule, length_list) : Appends A-A interactions for
 *                                                                          molecules of the same type to 
 *                                                                          the interaction list
 *                                                                          [A being dif atom types within
 *                                                                          the molecule] returns number
 *                                                                          added to the list of pair ints
 * ***************************************************************************************************/
int list_self_interactions(tW_word ** ilist, tW_molecule * mol, int* len_ilist)
{
  int length_added = 0;
  for (int i = 0; i < mol->n_apm; i++)
  {
    if (check_pairlist(mol->atom_types[i], mol->atom_types[i], ilist, *len_ilist) == FALSE)
    {
      strcpy(ilist[*len_ilist-1][0], mol->atom_types[i]);
      strcpy(ilist[*len_ilist-1][1], mol->atom_types[i]);

      *len_ilist= *len_ilist + 1;
      ilist = (tW_word **) erealloc(ilist, *len_ilist * sizeof(tW_word));
      ilist[*len_ilist-1] = (tW_word*) erealloc(ilist[*len_ilist-1], 2 * sizeof(tW_word));
      length_added++;
    }
  }
  return length_added;
}

/***********************************************************************************************
 * check_pairlist(word1, word2, pair_list, len_pair_list) : checks if the pair of word1 and word2
 *							    are in pair_list, returns true if so
 ***********************************************************************************************/
bool check_pairlist(tW_word typeA, tW_word typeB, tW_word **ilist, int len_list)
{
  bool found = FALSE;
  for (int i = 0; i < len_list; i++)
  {
   found = (
	     (strcmp(typeA, ilist[i][0]) == 0 && strcmp(typeB, ilist[i][1]) == 0 ) //If I have to do this again 
             ||								           //probably should just make a function
             (strcmp(typeA, ilist[i][1]) == 0 && strcmp(typeB, ilist[i][0]) == 0 ) 
	   );
   if (found == TRUE)
   {
     break;
   }
  }
  return found;
}

/***********************************************************************************************
 * check_list(int_list, list_len, my_int) : Checks int_list for my_int returns true if found
 * *******************************************************************************************/
bool check_list(int* list_of_int, int len, int my_int)
{
  bool found = FALSE;
  for (int i = 0; i < len; i++)
  {
    found = (list_of_int[i] == my_int);
    if (found == TRUE)
    {
      break;
    }
  }
return found;
}

/**********************************************************************************************************************
 * add_to_nb_list(interaction_list, molecule_to_check, current_list_length)  : Searches atom types within a molecule
 *                                                                             and adds pairs to running list of 
 *                                                                             interactions if permitted by exclusions
 *                                                                             updates list length accordingly
 *                                                                             returns # of interactions added to list
 *********************************************************************************************************************/
int add_to_nb_list(tW_word ** ilist, tW_molecule * mol, int* len_ilist)
{
  int count_added = 0;
  
  // look at all the atoms/sites in the molecule
  for (int i = 0; i < mol->n_apm; i++) // sit on an atom
  {
   for (int j = mol->n_apm-1; j>i; j--) // look at all the other atoms in this molecule we havent sat on yet 
   {
     if (check_list(mol->excls[i],mol->n_epa[i],j) == FALSE) // If we need to consider this atom pair
     {
       // check if we have the types in ilist alread
       if (check_pairlist(mol->atom_types[i], mol->atom_types[j], ilist, *len_ilist) == FALSE)
       {
         // We need to add these atom_types to the list and make room for more
         strcpy(ilist[*len_ilist-1][0], mol->atom_types[i]);
         strcpy(ilist[*len_ilist-1][1], mol->atom_types[j]);
         
	 *len_ilist= *len_ilist + 1;
         ilist = (tW_word **) erealloc(ilist, *len_ilist * sizeof(tW_word));
	 ilist[*len_ilist-1] = (tW_word*) erealloc(ilist[*len_ilist-1], 2 * sizeof(tW_word));
         count_added++;
       }  
     }
   }
  }
  return count_added;
}

/**************************************************************************************************************************************
 *sort_moltypes(system_top, singlemol_instances, multi_mol_instances, count_of_each) : lists moltypes with multiple occurances [many SOL]
 *                                                                                     and single instances [1 Protein]
 *                                                                                     in two different lists, 
 *                                                                                     tracks length of lists+1 in the dim_holder passed
 *************************************************************************************************************************************/
void sort_moltypes(tW_gmx_topology *top, int* one_list, int* multi_list, int* dim_holder)
{
 int multi_list_dim = 1;
 int one_list_dim = 1; 
 for (int i = 0; i < top->contents->mols.nr; ++i)
   {
     if (top->molecules[i].n_mols == 1)
     {
       //add it to the one_list and make room for more entries
       one_list[one_list_dim-1] = i;
       one_list_dim++;
       one_list = (int *) erealloc(one_list, one_list_dim * sizeof(int));
     }
     else 
     {
       //add it to the multi_list and make room for more entries
       multi_list[multi_list_dim-1] = i;
       multi_list_dim++;
       multi_list = (int *) erealloc(multi_list, multi_list_dim * sizeof(int));
     }
   }
 
 dim_holder[0] = one_list_dim;
 dim_holder[1] = multi_list_dim;
}

/*****************************************************************************************
 * brutal_count_nb_interactions(system_top) : returns the number of required nb interactions
 * 					      through brute force enumeration (sorry)
**************************************************************************************** */
int brutal_count_nb_interactions(tW_gmx_topology *top)
{
  int total = 0;
  
  // Lets seperate our topology into the single instance molecules and the duplicates
  int* lone_indices = (int*) emalloc(sizeof(int));
  int* multi_indices = (int*) emalloc(sizeof(int));
  int list_lengths[2];
  sort_moltypes(top, lone_indices, multi_indices, list_lengths);
  
  // Now we'll manually list the pairs and add to the list as we find unique pairs 
  tW_word ** inter_list = ecalloc(2000, sizeof(tW_word*)); // I am hard clearing alot of memory and hoping we don't have to go too much further
  inter_list[0] = ecalloc(2, sizeof(tW_word));
  int leni_list = 1;
  
  // count site interactions from long molecules that occur 1x 
  for (int i = 0; i < (list_lengths[0]-1); i++)
  {
   total = total + add_to_nb_list(inter_list, &(top->molecules[lone_indices[i]]), &leni_list);
  }
  
  // count site interactions from molecules that see another copy of themselves
  for (int j = 0; j < (list_lengths[1]-1); j++)
  {
    total = total + list_self_interactions(inter_list, &(top->molecules[multi_indices[j]]), &leni_list);
    total = total + list_cross_interactions(inter_list, &(top->molecules[multi_indices[j]]), &leni_list);
  }
  
  // Count Cross Molecular Interactions 
  total = total + list_intermol_interactions(inter_list, top, &leni_list); 
  return total;
}
/************************************************************************************
 * naive_exclusion_check(system_top) : checks if n_molecules > 1 for each moltype 
 *                                     in the system, returns true if there are
**************************************************************************************/
bool naive_exclusion_check(tW_gmx_topology *top)
{
   bool count_naiveway;
   int sum = 0;
   
   for (int i = 0; i < top->contents->mols.nr; ++i)
   {
     if (top->molecules[i].n_mols == 1)
     {
       sum++; // then we need to consider exclusions
     }
   }

   count_naiveway = (sum == 0);
   return count_naiveway;
}

/******************************************************************************************** 
 * count_required_nb(system_top) : determines the total nb interaction types needed based on top
 *                                 via either brute force listing or combinatorics
 *********************************************************************************************/
int count_required_nb(tW_gmx_topology *top)
{
  int total_nb = 0;
  tW_word** unique_interaction_list = ecalloc(1, sizeof(tW_word *)); // We won't use this in the function here, just need something to pass in list_nb_interactions
  // Count the easy case where we don't really care about bonded exclusions
  bool multiples_of_all_molecules = naive_exclusion_check(top); // We don't care about exclusions when we have multiple molecules of the same type
  if (multiples_of_all_molecules == TRUE)
  {
    //printf("Hooray, regular combinatorics work\n");
    total_nb = choose(top->n_atomtypes, 2) + top->n_atomtypes;
  }
  else 
  {
    //printf("This calls for something more complicated and PAINFUL\n");
    // If I had thought about this in advance, I would've counted and listed interactions at the same time
    total_nb = brutal_count_nb_interactions(top);
  }
  return total_nb;
}

/******************************************************************************************************
 * search_dihedrals(molecule, were_equiv_dihs_found, index_to_check) : loops through all ijkl dihs
 *                                                                   sets beq to true where equivalent 
 *                                                                   types found
 ******************************************************************************************************/
int search_dihedrals(tW_molecule * mol, bool * beq_found, int check)
{
  //Find the type of atoms involved in this dihedral
  tW_word nameA, nameB, nameC, nameD;
  strcpy(nameA,  mol->atom_types[mol->dih_ijkl[check][0]]);
  strcpy(nameB, mol->atom_types[mol->dih_ijkl[check][1]]); // Central
  strcpy(nameC, mol->atom_types[mol->dih_ijkl[check][2]]); // Pair Types
  strcpy(nameD, mol->atom_types[mol->dih_ijkl[check][3]]);
  //See if any other dihedral instances involve these same types
  for (int i = 0; i < mol->n_dihs; i++)
  {
  /*  if (
	( strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 && strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 ) ||
	( strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 && strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 )
       )
       {
    	 if (
            (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][0]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][3]]) == 0) ||
            (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][3]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][0]]) == 0)
            )
            {
              beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
            }
       }  */
  if ( (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][0]]) == 0 && strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 && 
        strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][3]]) == 0) ||
        (strcmp(nameA, mol->atom_types[mol->dih_ijkl[i][3]]) == 0 && strcmp(nameB, mol->atom_types[mol->dih_ijkl[i][2]]) == 0 &&
        strcmp(nameC, mol->atom_types[mol->dih_ijkl[i][1]]) == 0 && strcmp(nameD, mol->atom_types[mol->dih_ijkl[i][0]]) == 0)
     )
      {
        beq_found[i]=TRUE;
      }
  }
  return 1;
}

/**********************************************************************************************
 * get_n_distinct_dihedrals(system_top) : Counts the distinct dihedral potentials we need
 *                                        based on atomtype connectivity
 **********************************************************************************************/
int get_n_distinct_dihedrals(tW_gmx_topology *top)
{
  int total_dihedrals = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the dihedral before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_dihs, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_dihs; k++){ repeat_found[k] = FALSE;}

    //Search through dihedrals, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_dihs; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        total_dihedrals = total_dihedrals + search_dihedrals(&(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
      }
    }
  }
  return total_dihedrals;
}

/******************************************************************************************************
 * search_angles(molecule, were_equiv_angles_found, index_to_check) : loops through all ijk angles
 *                                                                    sets beq to true where equivalent 
 *                                                                    types found
 ******************************************************************************************************/
int search_angles(tW_molecule * mol, bool * beq_found, int check)
{
  //Find the type of atoms involved in this angle
  tW_word nameA, nameB, nameC;
  strcpy(nameA,  mol->atom_types[mol->angle_ijk[check][0]]);
  strcpy(nameB, mol->atom_types[mol->angle_ijk[check][1]]); // Central Atom
  strcpy(nameC, mol->atom_types[mol->angle_ijk[check][2]]);
  //See if any other instance of same angle type is present
  for (int i = 0; i < mol->n_angles; i++)
  {
    if (strcmp(nameB, mol->atom_types[mol->angle_ijk[i][1]]) == 0)
    {
      if (
	  (strcmp(nameA, mol->atom_types[mol->angle_ijk[i][0]]) == 0 && strcmp(nameC, mol->atom_types[mol->angle_ijk[i][2]]) == 0) ||
          (strcmp(nameA, mol->atom_types[mol->angle_ijk[i][2]]) == 0 && strcmp(nameC, mol->atom_types[mol->angle_ijk[i][0]]) == 0)
         )
         {
           beq_found[i]=TRUE;  // Mark all of the types that are duplicate as true
         }
    }
  }
  return 1;
}

/**********************************************************************************************
 * get_n_distinct_angles(system_top) : Counts the distinct angle potentials we need
 *                                     based on atomtype connectivity
 **********************************************************************************************/
int get_n_distinct_angles(tW_gmx_topology *top)
{
  int total_angles = 0;
  for (int i = 0; i < top->n_molecule_types; i++)
  {
    //Store if we've found the set of molecules in the angle before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_angles, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_angles; k++){ repeat_found[k] = FALSE;}
    
    //Search through angles, but not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_angles; j++)
    {
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
        total_angles = total_angles + search_angles(&(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
      }
    }
  }  
  return total_angles;
}

/******************************************************************************************
 * search_bonds(molecule, did_we_find_eq_set, index_to_check) : loops through all ij bonds
 *                                                              sets beq_found to true where
 *                                                              equivalent types found
 ******************************************************************************************/
int search_bonds(tW_molecule * mol, bool * beq_found, int check)
{
  //Find names we need to check
  char nameA[50], nameB[50];
  strcpy(nameA,  mol->atom_types[mol->bond_ij[check][0]]);
  strcpy(nameB, mol->atom_types[mol->bond_ij[check][1]]);
  
  //See if any other bonds involve these same two types
  for (int i = 0; i < mol->n_bonds; ++i)
  {
    if ( 
       (strcmp(nameA, mol->atom_types[mol->bond_ij[i][0]]) == 0 && strcmp(nameB, mol->atom_types[mol->bond_ij[i][1]]) == 0 ) 
       || 
       (strcmp(nameA, mol->atom_types[mol->bond_ij[i][1]]) == 0 && strcmp(nameB, mol->atom_types[mol->bond_ij[i][0]]) == 0)
       ){
	  beq_found[i]=TRUE; 
        }  // Mark all of the types that are duplicate as true
  }
  return 1;
}

/*********************************************************************************************
 * get_n_distinct_bonds(system_top) : Counts the distinct bond types needed based 
 * 				      on atomtype connectivity 
 *******************************************************************************************/
int get_n_distinct_bonds(tW_gmx_topology *top)
{
  int total_bonds = 0; 

  for (int i = 0; i < top->n_molecule_types; i++) // check bond total for each molecule type
  {
    // Store if we've checked a pair in the molecule before
    bool *repeat_found = (bool *) ecalloc(top->molecules[i].n_bonds, sizeof(bool));
    for (int k = 0; k < top->molecules[i].n_bonds; k++){ repeat_found[k] = FALSE;}
    
    // Search pairs, not the ones we already know have repeats
    for (int j = 0; j < top->molecules[i].n_bonds; j++)
    { 
      if (repeat_found[j] == FALSE) // count only attempts we havent found duplicates for yet
      {
	total_bonds = total_bonds + search_bonds(&(top->molecules[i]), repeat_found, j); // returns 1 when called, sets repeat_found elements to true
      } 
    }
  }
  return total_bonds; 
}

/*********************************************************************************************
 * choose(n,k) :  n choose k [as in combination math] (done recursively)
 ********************************************************************************************/
int choose(int n, int k)
{
  if (k == 0) return 1;
  return (n * choose(n-1,k-1))/k;
}

/*********************************************************************************************
 * calc_total_interactions(system_top) : calculates total number of interactions needed to
 *					 simulate topology from .btp file, this function is 
 *					 an alternative to contents->idef.ntypes since 
 *                                       that count includes AB AND BA interactions + simplifies 
 *                                       bead interactions if the ff that generated the .tpr file
 *                                       had multiple types with the identical functions of r
 **********************************************************************************************/
int calc_total_interactions(tW_gmx_topology *top)
{
  int total;
  // Count nb AB and AA ints - which is more complicated than I'd like 
  total = count_required_nb(top);
  printf("NONBONDED %d\n", total);
  
  // Count distinct bonded interactions
  total = total + get_n_distinct_bonds(top);
  printf("BONDED %d\n",get_n_distinct_bonds(top));
  
  // Count distinct angle interactions 
  total = total + get_n_distinct_angles(top);
  printf("ANGLES %d\n", get_n_distinct_angles(top));
  
  // Count distinct dihedral interactions
  total = total + get_n_distinct_dihedrals(top);
  printf("DIHEDRALS %d \n", get_n_distinct_dihedrals(top));

  total = total + get_n_impropernb(top);
  printf("LJ14 %d\n", get_n_impropernb(top));
  return total;
}

/*********************************************************************************************
 * write_par(filename, system_topology) : Prints information from the previously read in topology 
 * 					  in the format of a "par.txt" file. Output is intended 
 * 					  as template to avoid by hand interaction listing
 * 					  for simulation specific information I have left an 
 * 					  XXX for the user to fill in manually
 *********************************************************************************************/
void write_par(tW_word fnm, tW_gmx_topology *top)
{
  FILE *fp = open_file(fnm, 'w');
  //Header Info for user to fill
  fprintf(fp,"!# input file \n");
  fprintf(fp,"\n");
  fprintf(fp,"[Mode] GROMACS\n");
  fprintf(fp,"\n");
  fprintf(fp,"[Temperature]  XXX\n");
  fprintf(fp,"\n");
  fprintf(fp,"[Structures]    1\n");
  fprintf(fp,"inp.txt\n");
  fprintf(fp,"[End Structures]\n");

  // System Specific Info 
  fprintf(fp,"[Site_Types]    %d \n", top->n_atomtypes);
  for (int i = 0; i < top->n_atomtypes; i++)
  {
      fprintf(fp,"%8s\n", top->atom_type_names[i]);
  }
  fprintf(fp,"[End Site_Types]\n");

  //int total_inter_types = calc_total_interactions(top); Lets cut down on the number of times we have to loop through the top
  int total_inter_types, num_nbi, num_boni, num_angi, num_dihi, num_imnbi;
  
  num_nbi = count_required_nb(top);
  num_boni = get_n_distinct_bonds(top);
  num_angi = get_n_distinct_angles(top);
  num_dihi = get_n_distinct_dihedrals(top);
  num_imnbi = get_n_impropernb(top);

  total_inter_types = num_nbi + num_boni + num_angi + num_dihi + num_imnbi; 
  
  fprintf(fp,"[Inter_Types]   %d\n", total_inter_types);

  
  // Now the fun part, we need to list the appropriate types, luckily we know how many of each
  tW_word** pair_interaction_list = ecalloc(num_nbi, sizeof(tW_word*));
  tW_word** bonded_interaction_list = ecalloc(num_boni, sizeof(tW_word*));
  tW_word** angle_interaction_list = ecalloc(num_angi, sizeof(tW_word*));
  tW_word** dihedral_interaction_list = ecalloc(num_dihi, sizeof(tW_word*));
  tW_word** impropernb_interaction_list = ecalloc(num_imnbi, sizeof(tW_word*));


  for (int i = 0; i < num_nbi; i++) { pair_interaction_list[i] = ecalloc(2, sizeof(tW_word));}
  for (int i = 0; i < num_boni; i++) { bonded_interaction_list[i] = ecalloc(2, sizeof(tW_word));}
  for (int i = 0; i < num_angi; i++) { angle_interaction_list[i] = ecalloc(3, sizeof(tW_word));}
  for (int i = 0; i < num_dihi; i++) { dihedral_interaction_list[i] = ecalloc(4, sizeof(tW_word));}
  for (int i = 0; i < num_imnbi; i++) { impropernb_interaction_list[i] = ecalloc(2, sizeof(tW_word));}
  
  printf("Filling nb list, searching for %d interactions\n", num_nbi);
  fill_nb_list(top, pair_interaction_list, num_nbi);
  printf("Filling bonded list, searching for %d interactions\n", num_boni);
  fill_bond_list(top, bonded_interaction_list, num_boni);
  printf("Filling angle list, searching for %d interactions\n", num_angi);
  fill_angle_list(top, angle_interaction_list, num_angi);
  printf("Filling dihedral list, searching for %d interactions \n", num_dihi);
  fill_dihedral_list(top, dihedral_interaction_list, num_dihi);
  printf("Filling improper list, searching for %d interactions \n", num_imnbi);
  fill_impropernb_list(top, impropernb_interaction_list, num_imnbi);
  
  write_total_list(fp, pair_interaction_list, num_nbi, mNB);
  write_total_list(fp, bonded_interaction_list, num_boni, mBOND);
  write_total_list(fp, angle_interaction_list, num_angi, mANGLE);
  write_total_list(fp, dihedral_interaction_list, num_dihi, mDIH);
  write_total_list(fp, impropernb_interaction_list, num_imnbi, mIMNB);
  fprintf(fp, "[End Inter_Types]\n");

  write_interactions(fp, pair_interaction_list, num_nbi, mNB);
  write_interactions(fp, bonded_interaction_list, num_boni, mBOND);
  write_interactions(fp, angle_interaction_list, num_angi, mANGLE);
  write_interactions(fp, dihedral_interaction_list, num_dihi, mDIH);
  write_interactions(fp, impropernb_interaction_list, num_imnbi, mIMNB);
  


  fprintf(fp, "\n[TPR] 1\nYOUR_BTP_FILE\n [End TPR]");
  fclose(fp);
}
