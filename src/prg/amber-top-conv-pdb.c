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
#include <mol/molec.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* convert amber topology and coordinates to pdb molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates */
struct pdb_mol *amber_top_pdb(struct amber_top *t, double *c) {
  unsigned i,j,k,a0,a1,id=1;
  struct pdb_mol *p;
  /* pdb structure */
  p = mol_pdb_new();
  p->title = str_copy_new("The whole structure");
  p->remark = str_copy_new("Converted from Amber topology");
  /* residuii */
  p->n_res = t->pointers[AMBER_POINTER_NRES];
  p->res = mol_pdb_res_new(p->n_res);
  for (i=0; i<p->n_res; i++) {
    p->res[i].name = str_copy_new(t->res_names[i]);
    p->res[i].ter = amber_res_is_ter(t,i);
    p->res[i].id = i+1;
    /* atoms */
    amber_res_atom_id(t,i,&a0,&a1);
    p->res[i].n_atoms = a1-a0;
    p->res[i].atom = mol_pdb_atom_new(p->res[i].n_atoms);
    for (j=0; j<p->res[i].n_atoms; j++) {
      for (k=0; k<3; k++)
        p->res[i].atom[j].coord[k] = c[3*(a0+j)+k];
      p->res[i].atom[j].name = str_copy_new(t->atom_names[a0+j]);
      p->res[i].atom[j].charge = t->charge[a0+j]/AMBER_TOP_CHARGE;
      p->res[i].atom[j].id = id++;
      }
    }
  return(p);
  }

/* convert amber topology and coordinates to pdb molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates
   r - residuum ID */
struct pdb_mol *amber_top_pdb_r(struct amber_top *t, double *c, unsigned r) {
  unsigned j,k,a0,a1,id = 1;
  char text[1024];
  struct pdb_mol *p = NULL;
  if (!r)
    p = amber_top_pdb(t,c);
  else {
    p = mol_pdb_new();
    sprintf(text,"Residuum #%d: %s",r,t->res_names[r-1]);
    p->title = str_copy_new(text);
    p->remark = str_copy_new("Converted from Amber topology");
    p->n_res = 1;
    p->res = mol_pdb_res_new(p->n_res);
    p->res[0].name = str_copy_new(t->res_names[r-1]);
    p->res[0].ter = amber_res_is_ter(t,r-1);
    p->res[0].id = r;
    amber_res_atom_id(t,r-1,&a0,&a1);
    p->res[0].n_atoms = a1-a0;
    p->res[0].atom = mol_pdb_atom_new(p->res[0].n_atoms);
    for (j=0; j<p->res[0].n_atoms; j++) {
      for (k=0; k<3; k++)
        p->res[0].atom[j].coord[k] = c[3*(a0+j)+k];
      p->res[0].atom[j].name = str_copy_new(t->atom_names[a0+j]);
      p->res[0].atom[j].charge = t->charge[a0+j]/AMBER_TOP_CHARGE;
      p->res[0].atom[j].id = id++;
      }
    }
  return(p);
  }

/* convert amber topology and coordinates to pdb molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates
   v - array with residuum IDs
   n - number of residui */
struct pdb_mol *amber_top_pdb_rv(struct amber_top *t, double *c,
  unsigned *v, unsigned n) {
  unsigned i,j,k,a0,a1,id = 1;
  char range[1024],text[2048];
  struct pdb_mol *p = NULL;
  if (!v || !n)
    p = amber_top_pdb(t,c);
  else {
    p = mol_pdb_new();
    str_unum_range(range,v,n);
    sprintf(text,"Residui #%s",range);
    p->title = str_copy_new(text);
    p->remark = str_copy_new("Converted from Amber topology");
    p->n_res = n;
    p->res = mol_pdb_res_new(p->n_res);
    for (i=0; i<n; i++) {
      p->res[i].name = str_copy_new(t->res_names[v[i]-1]);
      p->res[i].ter = amber_res_is_ter(t,v[i]-1);
      p->res[i].id = v[i];
      amber_res_atom_id(t,v[i]-1,&a0,&a1);
      p->res[i].n_atoms = a1-a0;
      p->res[i].atom = mol_pdb_atom_new(p->res[i].n_atoms);
      for (j=0; j<p->res[i].n_atoms; j++) {
        for (k=0; k<3; k++)
          p->res[i].atom[j].coord[k] = c[3*(a0+j)+k];
        p->res[i].atom[j].name = str_copy_new(t->atom_names[a0+j]);
        p->res[i].atom[j].charge = t->charge[a0+j]/AMBER_TOP_CHARGE;
        p->res[i].atom[j].id = id++;
        }
      }
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
