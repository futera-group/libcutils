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
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* return internal code of orbital type (s,p,d,f)
 
   s - orbital label */
short cpmd_bs_id(char *s) {
  if (str_compare(s,"S"))
    return(CPMD_BASIS_S);
  if (str_compare(s,"P"))
    return(CPMD_BASIS_P);
  if (str_compare(s,"D"))
    return(CPMD_BASIS_D);
  if (str_compare(s,"F"))
    return(CPMD_BASIS_F);
  return(CPMD_BASIS_S);
  }

/* convert inter code of orbital type to its symbol (s,p,d,f)
   
   s - orbital label
   n - orbital id
   b - subtype id */
char* cpmd_bs_sym(short n, unsigned b) {
  static char s[10]="\0";
  char sym[10];
  switch (n) {
    case CPMD_BASIS_S: sprintf(s,"S"); break;
    case CPMD_BASIS_P: sprintf(s,"P"); break;
    case CPMD_BASIS_D: sprintf(s,"D"); break;
    case CPMD_BASIS_F: sprintf(s,"F"); break;
    default: sprintf(s,"X");
    }
  if (n!=CPMD_BASIS_S) {
    sprintf(sym,"%s",s);
    switch (b) {
      case 0: sprintf(s,"%s%s",sym,"x"); break;
      case 1: sprintf(s,"%s%s",sym,"y"); break;
      case 2: sprintf(s,"%s%s",sym,"z"); break;
      }
    }
  return(s);
  }

/* create array with atom orbital types for specified type of atom
 
   d - pointer to cpmd data struct
   q - queue with atom orbital IDs
   t - ID of atomic type */
void cpmd_bs_set(struct cpmd_dat *d, struct queue *q, unsigned t) {
  unsigned id = 0,n_tot = 0;
  short n;
  if (d->type[t].ao_id)
    d->type[t].ao_id = vec_sifree(d->type[t].ao_id);
  d->type[t].n_ao_ids = q->num;
  d->type[t].ao_id = vec_sialloc(d->type[t].n_ao_ids);
  while (q->num) {
    queue_siget(q,&n);
    d->type[t].ao_id[id] = n;
    n_tot += n;
    id++;
    }
  d->type[t].n_orbs = n_tot;
  }

/* -------------------------------------------------------------------------- */
