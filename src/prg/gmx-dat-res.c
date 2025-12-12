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

#include <stdio.h>
#include <stdlib.h>
#include <cmn/list.h>
#include <cmn/types.h>
#include <cmn/vector.h>
#include <mol/distance.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Calculate COM of given residuum

   g  - gromacs data structure
   r  - residuum data structure
   a0 - atomic offset (ID of the first atom of the residuum)
   c  - array for storing the COM coordinates */
void gmx_dat_res_com(struct gmx_dat *g, struct gmx_res *r, unsigned a0, 
  double *c) {
  unsigned i,j;
  double m_tot = 0.0;
  vec_fset(c,0.0,3);
  for (i=0; i<r->n_atoms; i++) {
    for (j=0; j<3; j++)
      c[j] += r->atom[i].mass*g->crd[a0+i][j];
    m_tot += r->atom[i].mass;
    }
  for (i=0; i<3; i++)
    c[i] /= m_tot;
  }

/* Return IDs of residui within given distance of specific residui

   g - gromacs data structure
   v - array with residuum IDs
   n - number of residui */
unsigned* gmx_dat_res_get_list_dr(struct gmx_dat *g, unsigned *v, unsigned n,
  double d, unsigned *nr) {
  unsigned i,j,k,l,ia,ir,ii,*rr;
  int id;
  double r,com[3];
  struct list *t,*w;
  struct ldata *p;
  /* list of specified residuum atoms */
  t = list_alloc();
  for (i=0,ia=0,ir=0; i<g->n_frags; i++)
    for (j=0; j<g->frag[i].n_rep; j++)
      for (k=0; k<g->frag[i].mol->n_resids; k++,ir++) {
        id = vec_ufind(v,ir+1,n);
        for (l=0; l<g->frag[i].mol->res[k].n_atoms; l++,ia++)
          if (id >= 0)
            list_add_end(t,&ia,TYPE_UINT);
        }
  /* list of residues within the specified distance */
  w = list_alloc();
  for (i=0,ia=0,ir=0; i<g->n_frags; i++)
    for (j=0; j<g->frag[i].n_rep; j++)
      for (k=0; k<g->frag[i].mol->n_resids; k++,ir++) {
        if (vec_ufind(v,ir+1,n) < 0) {
          /* center of mass */
          gmx_dat_res_com(g,g->frag[i].mol->res+k,ia,com);
          /* distance to specified atoms */
          for (p=t->first; p; p=p->l_next) {
            ii = *((unsigned*)(p->l_data));
            r = dist_r(com,g->crd[ii],3);
            if (r <= d) {
              list_add_end(w,&ir,TYPE_UINT);
              break;
              }
            }
          }
        ia += g->frag[i].mol->res[k].n_atoms;
        }
  list_free(t,free);
  /* convert list to array */
  *nr = w->num;
  rr = vec_ualloc(*nr);
  for (p=w->first,i=0; p; p=p->l_next,i++)
    rr[i] = *((unsigned*)(p->l_data));
  list_free(w,free);
  return rr;
  }

/* -------------------------------------------------------------------------- */
