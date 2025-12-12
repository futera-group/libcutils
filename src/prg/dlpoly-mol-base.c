/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                28/01/2017                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address:             University College London (UCL)                      * 
 *                       Department of Physics and Astronomy                  * 
 *                       Gower Street, London WC1E 6BT                        * 
 *                       United Kingdom                                       * 
 *                                                                            * 
 *  E-Mail:              z.futera@ucl.ac.uk                                   * 
 *                                                                            * 
\******************************************************************************/

#include <stdlib.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY molecular data
   
   n - number of the molecules */
struct dlpoly_mol *dlpoly_mol_new(unsigned n) {
  unsigned i;
  struct dlpoly_mol *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_mol));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY molecular data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].name = NULL;
    d[i].n_molecules = 0;
    d[i].n_atoms = 0;
    d[i].atom = NULL;
    d[i].n_shells = 0;
    d[i].shell = NULL;
    d[i].n_constraints = 0;
    d[i].cnst = NULL;
    d[i].n_pmfs = 0;
    d[i].pmf = NULL;
    d[i].n_rigids = 0;
    d[i].rig = NULL;
    d[i].n_teths = 0;
    d[i].teth = NULL;
    d[i].n_bonds = 0;
    d[i].bond = NULL;
    d[i].n_angles = 0;
    d[i].angle = NULL;
    d[i].n_dihedrals = 0;
    d[i].dihed = NULL;
    d[i].n_inversions = 0;
    d[i].inv = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY molecular data 
 
   d - the molecular data
   n - number of molecules */
struct dlpoly_mol *dlpoly_mol_free(struct dlpoly_mol *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].name = str_free(d[i].name);
      d[i].atom = dlpoly_atom_free(d[i].atom,d[i].n_atoms);
      d[i].shell = dlpoly_prm_shell_free(d[i].shell);
      d[i].cnst = dlpoly_prm_cnst_free(d[i].cnst);
      d[i].pmf = dlpoly_prm_pmf_free(d[i].pmf);
      d[i].rig = dlpoly_prm_rig_free(d[i].rig,d[i].n_rigids);
      d[i].teth = dlpoly_prm_teth_free(d[i].teth);
      d[i].bond = dlpoly_prm_bond_free(d[i].bond);
      d[i].angle = dlpoly_prm_angle_free(d[i].angle);
      d[i].dihed = dlpoly_prm_dihed_free(d[i].dihed);
      d[i].inv = dlpoly_prm_inv_free(d[i].inv);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
