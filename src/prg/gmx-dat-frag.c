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

#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Add molecular coordinates to the system
 
   d  - gromacs data structure 
   r  - array of fragment data structures
   nr - number of molecules
   c  - array of coordinates
   nc - number of coordinates */
void gmx_dat_frag_add(struct gmx_dat *d, struct gmx_frag *r, unsigned nr, 
  double **c, unsigned nc) {
  unsigned i,nt;
  double **x;
  struct gmx_frag *t;
  /* new array of fragments */
  nt = d->n_frags + nr;
  t = gmx_frag_new(nt);
  for (i=0; i<d->n_frags; i++)
    gmx_frag_copy(t+i,d->frag+i);
  for (i=0; i<nr; i++) {
    gmx_frag_copy(t+d->n_frags+i,r+i);
    t[d->n_frags+i].mol = gmx_dat_mol_get(d,r[i].name);
    }
  /* replace the array */
  d->frag = gmx_frag_free(d->frag,d->n_frags);
  d->n_frags = nt;
  d->frag = t;
  /* new array of coordinates */
  nt = d->n_atoms + nc;
  x = mat_falloc(nt,3);
  for (i=0; i<d->n_atoms; i++)
    vec_fcopy(x[i],d->crd[i],3);
  for (i=0; i<nc; i++)
    vec_fcopy(x[d->n_atoms+i],c[i],3);
  /* replace the array */
  d->crd = mat_ffree(d->crd,d->n_atoms);
  d->crd = x;
  /* update number of atoms/residues */
  gmx_dat_update_nums(d);
  if (d->n_atoms!=nt)
    msg_error_f("inconsistent number of atoms (%d) and coordinates (%d)",
      1,d->n_atoms,nt);
  }

/* Delete specified molecular fragments from the system 
 
   d  - gromacs data structure
   id - array of frament IDs
   n  - number of fragments to delete */
void gmx_dat_frag_del(struct gmx_dat *d, unsigned *id, unsigned n) {
  unsigned i,j,ia,it,ix,nt,nc,nn;
  short *a;
  double **x;
  struct gmx_frag *t;
  /* array of delete/stay indicators */
  a = vec_sialloc(d->n_frags);
  vec_siset(a,1,d->n_frags);
  for (i=0; i<n; i++)
    a[id[i]] = 0;
  /* number of fragments/coordinates after deletion */
  for (i=0,nc=0,nt=0; i<d->n_frags; i++)
    if (a[i]) {
      nc += d->frag[i].n_rep*d->frag[i].mol->n_atoms;
      nt++;
      }
  /* delete all */
  if (nt==0) {
    d->frag = gmx_frag_free(d->frag,d->n_frags);
    d->n_frags = 0;
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->n_atoms = 0;
    d->n_resids = 0;
    }
  /* delete selected */
  else if (nt<d->n_frags) {
    t = gmx_frag_new(nt);
    x = mat_falloc(nc,3);
    for (i=0,ia=0,it=0,ix=0; i<d->n_frags; i++) {
      nn = d->frag[i].n_rep*d->frag[i].mol->n_atoms;
      if (a[i]) {
        for (j=0; j<nn; j++)
          vec_fcopy(x[ix++],d->crd[ia++],3);
        gmx_frag_copy(t+it++,d->frag+i);
        }
      else
        ia += nn;
      }
    d->frag = gmx_frag_free(d->frag,d->n_frags);
    d->n_frags = nt;
    d->frag = t;
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->crd = x;
    gmx_dat_update_nums(d);
    if (d->n_atoms!=nc)
      msg_error_f("inconsistent number of atoms (%d) and coordinates (%d)",
        1,d->n_atoms,nc);
    }
  /* clean memory */
  vec_sifree(a);
  }

/* Return fragment ID of specified atom in the system

   d - gromacs data structure
   ir - ID of residuum in molecule (optional)
   ia - ID of atom in residuum (optional)
   id - ID of the atom in system */
unsigned gmx_dat_frag_get_id(struct gmx_dat *d, unsigned *ir, unsigned *ia,
  unsigned id) {
  unsigned i,j,k,l,is;
  for (i=0,is=0; i<d->n_frags; i++)
    for (j=0; j<d->frag[i].n_rep; j++)
      for (k=0; k<d->frag[i].mol->n_resids; k++) 
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,is++) {
          if (is==id) {
            if (ir)
              (*ir) = k;
            if (ia)
              (*ia) = l;
            return(i);
            }
          }
  msg_error_f("atom ID %d not found in the system",1,id+1);
  return(0);
  }

/* -------------------------------------------------------------------------- */

/* Compact fragment array in gromacs data structure (use repetitions)
 
   d - gromacs data array */
void gmx_dat_frag_compact(struct gmx_dat *d) {
  unsigned i,j,n;
  struct queue *q;
  struct gmx_frag *x,*t;
  /* find duplicates */
  i = 0;
  q = queue_alloc();
  while (i<d->n_frags) {
    n = d->frag[i].n_rep;
    for (j=i+1; j<d->n_frags; j++) {
      if (!gmx_frag_compare(d->frag+i,d->frag+j))
        break;
      n += d->frag[j].n_rep;
      }
    x = gmx_frag_copy_new(d->frag+i,1);
    x->n_rep = n;
    queue_add(q,x);
    i = j;
    }
  /* create new array */
  n = q->num;
  t = gmx_frag_new(n);
  for (i=0; i<n; i++) {
    x = queue_get(q);
    gmx_frag_copy(t+i,x);
    gmx_frag_free(x,1);
    }
  /* replace the array */
  d->frag = gmx_frag_free(d->frag,d->n_frags);
  d->n_frags = n;
  d->frag = t;
  /* clean memory */
  queue_free(q);
  }

/* Expand fragment array in gromacs data structure (no repetitions)
 
   d - gromacs data array */
void gmx_dat_frag_expand(struct gmx_dat *d) {
  unsigned i,j,n;
  short found;
  struct queue *q;
  struct gmx_frag *x,*t;
  /* check whether there are any repetitions */
  for (i=0,found=0; i<d->n_frags; i++)
    if (d->frag[i].n_rep>1) {
      found = 1;
      break;
      }
  if (!found)
    return;
  /* expand duplicates */
  q = queue_alloc();
  for (i=0; i<d->n_frags; i++)
    for (j=0; j<d->frag[i].n_rep; j++) {
      x = gmx_frag_copy_new(d->frag+i,1);
      x->n_rep = 1;
      queue_add(q,x);
      }
  /* create new array */
  n = q->num;
  t = gmx_frag_new(n);
  for (i=0; i<n; i++) {
    x = queue_get(q);
    gmx_frag_copy(t+i,x);
    gmx_frag_free(x,1);
    }
  /* replace the array */
  d->frag = gmx_frag_free(d->frag,d->n_frags);
  d->n_frags = n;
  d->frag = t;
  /* clean memory */
  queue_free(q);
  }

/* -------------------------------------------------------------------------- */
