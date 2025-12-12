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

#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Set residuum IDs and sort bonding terms to inter/intra residual ones
  
   d   - molecular data structure */
void gmx_mol_set_bonding(struct gmx_mol *m) {
  unsigned i,j,k,id,*n,*nr,nt,n_terms = 10;
  short term[10] = {GMX_BTERM_BOND,GMX_BTERM_ANGL,GMX_BTERM_DIHE,
    GMX_BTERM_IMPR,GMX_BTERM_PAIR,GMX_BTERM_NBPR,GMX_BTERM_CNST,
    GMX_BTERM_STTL,GMX_BTERM_CMAP,GMX_BTERM_EXCL};
  short inter;
  struct gmx_bond **b,**br,*bt;
  struct queue *qm,**qr;
  /* bonding term */
  for (i=0; i<n_terms; i++) {
    b = gmx_top_bond_get(m->top,term[i],&n);
    if ((*b) && (*n)) {
      /* storage for term IDs */
      qm = queue_alloc();
      qr = vec_talloc(sizeof(struct queue*),m->n_resids);
      for (j=0; j<m->n_resids; j++)
        qr[j] = queue_alloc();
      for (j=0; j<(*n); j++) {
        /* fix atom/residuum IDs and atom names */
        (*b)[j].ir = vec_ufree((*b)[j].ir);
        (*b)[j].ir = vec_ualloc((*b)[j].n_atoms);
        (*b)[j].name = vec_talloc(sizeof(char*),(*b)[j].n_atoms);
        for (k=0; k<(*b)[j].n_atoms; k++) {
          (*b)[j].ir[k] = gmx_mol_res_get_id(m,&id,(*b)[j].ia[k]);
          (*b)[j].ia[k] = id;
          (*b)[j].name[k] = str_copy_new(m->res[(*b)[j].ir[k]].atom[id].name);
          }
        /* inter/intra-residual term sorting */
        for (k=1,inter=0; k<(*b)[j].n_atoms; k++)
          if ((*b)[j].ir[k]!=(*b)[j].ir[0]) {
            inter = 1;
            break;
            }
        if (inter)
          queue_uadd(qm,j);
        else
          queue_uadd(qr[(*b)[j].ir[0]],j);
        }
      /* intra-residual terms */
      for (j=0; j<m->n_resids; j++) {
        if (qr[j]->num) {
          br = gmx_top_bond_get(m->res[j].top,term[i],&nr);
          nt = (*nr)+qr[j]->num;
          bt = gmx_bond_new(nt);
          for (k=0; k<(*nr); k++)
            gmx_bond_copy(bt+k,(*br)+k);
          for (k=(*nr); k<nt; k++) {
            queue_uget(qr[j],&id);
            gmx_bond_copy(bt+k,(*b)+id);
            }
          (*br) = bt;
          (*nr) = nt;
          br = gmx_top_bond_get(m->res[j].top,term[i],&nr);
          }
        queue_free(qr[j]);
        }
      vec_tfree(qr);
      /* inter-residual terms */
      nt = qm->num;
      bt = gmx_bond_new(nt);
      for (j=0; j<nt; j++) {
        queue_uget(qm,&id);
        gmx_bond_copy(bt+j,(*b)+id);
        }
      (*b) = gmx_bond_free(*b,*n);
      (*n) = nt;
      (*b) = bt;
      queue_free(qm);
      }
    }
  }

/* -------------------------------------------------------------------------- */
