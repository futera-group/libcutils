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

#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/molec.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Convert gromacs data structure to XYZ molecular format

   d - gromacs data structure */
struct xyz_mol* gmx_dat_conv_to_xyz(struct gmx_dat *d) {
  unsigned i,j,k,l,ia;
  struct xyz_mol *x;
  /* structure */
  x = mol_xyz_new();
  x->title = str_copy_new(d->title);
  /* atoms */
  x->n_atoms = d->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  for (i=0,ia=0; i<d->n_frags; i++) {
    for (j=0; j<d->frag[i].n_rep; j++) {
      for (k=0; k<d->frag[i].mol->n_resids; k++) {
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
          x->atom[ia].name = str_copy_new(d->frag[i].mol->res[k].atom[l].name);
          x->atom[ia].num = d->frag[i].mol->res[k].atom[l].num;
          vec_fcopy(x->atom[ia].coord,d->crd[ia],3);
          }
        }
      }
    }
  return(x);
  }

/* Convert gromacs data structure to XYZ molecular format

   d - gromacs data structure
   v - array with residuum IDs
   n - number of residui */
struct xyz_mol* gmx_dat_conv_to_xyz_rv(struct gmx_dat *d, 
  unsigned *v, unsigned n) {
  unsigned i,j,k,l,ia,ir,id;
  char range[1024],text[2048];
  struct xyz_mol *x;
  if (!v || !n)
    x = gmx_dat_conv_to_xyz(d);
  else {
    x = mol_xyz_new();
    /* title */
    str_unum_range(range,v,n);
    sprintf(text,"Residui #%s",range);
    x->title = str_copy_new(text);
    /* number of atoms */
    x->n_atoms = 0;
    for (i=0,ir=0; i<d->n_frags; i++)
      for (j=0; j<d->frag[i].n_rep; j++)
        for (k=0; k<d->frag[i].mol->n_resids; k++,ir++)
          if (vec_ufind(v,ir+1,n) >= 0)
            x->n_atoms += d->frag[i].mol->res[k].n_atoms;
    /* atomic data */
    x->atom = mol_xyz_atom_new(x->n_atoms);
    for (i=0,ia=0,ir=0,id=0; i<d->n_frags; i++)
      for (j=0; j<d->frag[i].n_rep; j++) 
        for (k=0; k<d->frag[i].mol->n_resids; k++,ir++)
          for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) 
            if (vec_ufind(v,ir+1,n) >= 0) {
              x->atom[id].name =  
                str_copy_new(d->frag[i].mol->res[k].atom[l].name);
              x->atom[id].num = d->frag[i].mol->res[k].atom[l].num;
              vec_fcopy(x->atom[id].coord,d->crd[ia],3);
              id++;
              }
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
