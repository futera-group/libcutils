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

/* Convert gromacs data structure to PDB molecular format

   d - gromacs data structure */
struct pdb_mol* gmx_dat_conv_to_pdb(struct gmx_dat *d) {
  unsigned i,j,k,l,ia,ir;
  struct pdb_mol *p;
  /* structure */
  p = mol_pdb_new();
  p->title = str_copy_new(d->title);
  /* residues */
  p->n_res = d->n_resids; 
  p->res = mol_pdb_res_new(p->n_res);
  for (i=0,ia=0,ir=0; i<d->n_frags; i++) {
    for (j=0; j<d->frag[i].n_rep; j++) {
      for (k=0; k<d->frag[i].mol->n_resids; k++,ir++) {
        p->res[ir].name = str_copy_new(d->frag[i].mol->res[k].name);
        p->res[ir].ter = ((k+1)==d->frag[i].mol->n_resids ? 1 : 0);
        p->res[ir].id = ir+1;
        /* atoms */
        p->res[ir].n_atoms = d->frag[i].mol->res[k].n_atoms;
        p->res[ir].atom = mol_pdb_atom_new(p->res[ir].n_atoms);
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
          mol_pdb_atom_name_pos(p,ir,l,d->frag[i].mol->res[k].atom[l].name);
          p->res[ir].atom[l].id = ia+1;
          vec_fcopy(p->res[ir].atom[l].coord,d->crd[ia],3);
          p->res[ir].atom[l].charge = d->frag[i].mol->res[k].atom[l].charge;
          }
        }
      }
    }
  return(p);
  }

/* Convert gromacs data structure to PDB molecular format

   d - gromacs data structure
   v - array with residuum IDs
   n - number of residui */
struct pdb_mol* gmx_dat_conv_to_pdb_rv(struct gmx_dat *d, 
  unsigned *v, unsigned n) {
  unsigned i,j,k,l,ia,ir,id;
  char range[1024],text[2048];
  struct pdb_mol *p;
  if (!v || !n)
    p = gmx_dat_conv_to_pdb(d);
  else {
    p = mol_pdb_new();
    /* title */
    str_unum_range(range,v,n);
    sprintf(text,"Residui #%s",range);
    p->title = str_copy_new(text);
    /* residues */
    p->n_res = n;
    p->res = mol_pdb_res_new(p->n_res);
    for (i=0,ia=0,ir=0,id=0; i<d->n_frags; i++)
      for (j=0; j<d->frag[i].n_rep; j++)
        for (k=0; k<d->frag[i].mol->n_resids; k++,ir++) {
          if (vec_ufind(v,ir+1,n) >= 0) {
            p->res[id].name = str_copy_new(d->frag[i].mol->res[k].name);
            p->res[id].ter = ((k+1)==d->frag[i].mol->n_resids ? 1 : 0);
            p->res[id].id = ir+1;
            /* atoms */
            p->res[id].n_atoms = d->frag[i].mol->res[k].n_atoms;
            p->res[id].atom = mol_pdb_atom_new(p->res[id].n_atoms);
            for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
              mol_pdb_atom_name_pos(p,id,l,d->frag[i].mol->res[k].atom[l].name);
              p->res[id].atom[l].id = ia+1;
              vec_fcopy(p->res[id].atom[l].coord,d->crd[ia],3);
              p->res[id].atom[l].charge = d->frag[i].mol->res[k].atom[l].charge;
              }
            id++;
            }
          else
            ia += d->frag[i].mol->res[k].n_atoms;
          }
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
