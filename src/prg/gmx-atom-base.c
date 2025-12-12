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

/* Allocate memory for array of atoms

   n - number of atoms */
struct gmx_atom* gmx_atom_new(unsigned n) {
  unsigned i;
  struct gmx_atom *d = NULL;
  /* memory allocation */
  d = (struct gmx_atom*)malloc(n*sizeof(struct gmx_atom));
  if (!d)
    msg_error("cannot allocate memory for atomic data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].particle = 0;
    d[i].id = 0;
    d[i].num = 0;
    d[i].name = NULL;
    d[i].type = NULL;
    d[i].charge = 0.0;
    d[i].mass = 0.0;
    d[i].sigma = 0.0;
    d[i].epsilon = 0.0;
    }
  return(d);
  }

/* Clean memory allocated in atom structure

   d - the data structure */
void gmx_atom_clean(struct gmx_atom *d) {
  d->particle = 0;
  d->id = 0;
  d->num = 0;
  d->name = str_free(d->name);
  d->type = str_free(d->type);
  d->charge = 0.0;
  d->mass = 0.0;
  d->sigma = 0.0;
  d->epsilon = 0.0;
  }

/* Clean memory allocated for atom array

   d - the data array
   n - length of the array */
struct gmx_atom* gmx_atom_free(struct gmx_atom *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      gmx_atom_clean(d+i);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy atomic data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_atom_copy(struct gmx_atom *d1, struct gmx_atom *d2) {
  if (d1 && d2) {
    d1->particle = d2->particle;
    d1->id = d2->id;
    d1->num = d2->num;
    d1->name = str_copy_new(d2->name);
    d1->type = str_copy_new(d2->type);
    d1->charge = d2->charge;
    d1->mass = d2->mass;
    d1->sigma = d2->sigma;
    d1->epsilon = d2->epsilon;
    }
  }

/* Create copy of atom data structure array

   d - the data array
   n - length of the array */
struct gmx_atom* gmx_atom_copy_new(struct gmx_atom *d, unsigned n) {
  unsigned i;
  struct gmx_atom *t = NULL;
  if (d) {
    t = gmx_atom_new(n);
    for (i=0; i<n; i++)
      gmx_atom_copy(t+i,d+i);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* Compare two atom data structures
 
   a1,a2 - the data structures */
short gmx_atom_compare(struct gmx_atom *a1, struct gmx_atom *a2) {
  /* number */
  if (a1->num!=a2->num)
    return(0);
  /* particle */
  if (a1->particle!=a2->particle)
    return(0);
  /* name */
  if (!str_compare(a1->name,a2->name))
    return(0);
  /* type */
  if (!str_compare(a1->type,a2->type))
    return(0);
  /* charge */
  if (a1->charge!=a2->charge)
    return(0);
  /* mass */
  if (a1->mass!=a2->mass)
    return(0);
  return(1);
  }

/* -------------------------------------------------------------------------- */
