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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for array of fragments

   n - number of fragments */
struct gmx_frag* gmx_frag_new(unsigned n) {
  unsigned i;
  struct gmx_frag *d = NULL;
  /* memory allocation */
  d = (struct gmx_frag*)malloc(n*sizeof(struct gmx_frag));
  if (!d)
    msg_error("cannot allocate memory for fragment data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].type = 0;
    d[i].name = NULL;
    d[i].n_rep = 0;
    d[i].mol = NULL;
    }
  return(d);
  }

/* Clean memory allocated in fragment structure

   d - the data structure */
void gmx_frag_clean(struct gmx_frag *d) {
  d->name = str_free(d->name);
  d->n_rep = 0;
  d->mol = NULL;
  }

/* Clean memory allocated for fragment array

   d - the data array
   n - length of the array */
struct gmx_frag* gmx_frag_free(struct gmx_frag *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      gmx_frag_clean(d+i);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy fragment data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_frag_copy(struct gmx_frag *d1, struct gmx_frag *d2) {
  if (d1 && d2) {
    d1->type = d2->type;
    d1->name = str_copy_new(d2->name);
    d1->n_rep = d2->n_rep;
    d1->mol = d2->mol;
    }
  }

/* Create copy of fragment data structure array

   d - the data array
   n - length of the array */
struct gmx_frag* gmx_frag_copy_new(struct gmx_frag *d, unsigned n) {
  unsigned i;
  struct gmx_frag *t = NULL;
  if (d) {
    t = gmx_frag_new(n);
    for (i=0; i<n; i++)
      gmx_frag_copy(t+i,d+i);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* Compare two fragment data structures
 
   d1,d2 - the data structures */
short gmx_frag_compare(struct gmx_frag *d1, struct gmx_frag *d2) {
  /* name */
  if (!str_compare(d1->name,d2->name))
    return(0);
  /* pointer to molecule */
  if (d1->mol!=d2->mol)
    return(0);
  return(1);
  }

/* -------------------------------------------------------------------------- */
