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
#include <cmn/vector.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for basis set struct */
struct basis* basis_new(void) {
  struct basis *b = NULL;
  /* memory allocation */
  b = (struct basis*)malloc(sizeof(struct basis));
  if (!b)
    msg_error("cannot allocate basis set struct",1);
  /* initialization */
  b->n_bfce = 0;
  b->n_ibfce = 0;
  b->n_centers = 0;
  b->n_cont_shells = 0;
  b->n_prim_shells = 0;
  b->n_pure_d = 0;
  b->n_pure_f = 0;
  b->max_ang_mom = 0;
  b->max_bfce_cont = 0;
  b->center = NULL;
  return(b);
  }

/* allocate memory for basis set shell struct
 
   n - number of the shells */
struct basis_shell* basis_shell_new(unsigned n) {
  struct basis_shell *b = NULL;
  unsigned i;
  if (n) {
    /* memory allocation */
    b = (struct basis_shell*)malloc(n*sizeof(struct basis_shell));
    if (!b)
      msg_error("cannot allocate basis set shell struct",1);
    /* initialization */
    for (i=0; i<n; i++) {
      b[i].type = 0;
      b[i].bf = 0;
      b[i].n_bfce = 0;
      b[i].n_prim = 0;
      b[i].exp = NULL;
      b[i].cf1 = NULL;
      b[i].cf2 = NULL;
      }
    }
  return(b);
  }

/* allocate memory for basis set atom struct
 
   n - number of the atoms */
struct basis_center* basis_center_new(unsigned n) {
  struct basis_center *b = NULL;
  unsigned i,j;
  if (n) {
    /* memory allocation */
    b = (struct basis_center*)malloc(n*sizeof(struct basis_center));
    if (!b)
      msg_error("cannot allocate basis set center struct",1);
    /* initialization */
    for (i=0; i<n; i++) {
      b[i].type = 0;
      b[i].bf = 0;
      b[i].n_bfce = 0;
      b[i].n_shells = 0;
      b[i].shell = NULL;
      for (j=0; j<3; j++)
        b[i].coord[j] = 0.0;
      }
    }
  return(b);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for basis set
 
   b - pointer to basis set struct */
struct basis* basis_free(struct basis *b) {
  if (b) {
    b->center = basis_center_free(b->center,b->n_centers);
    free(b);
    }
  return(NULL);
  }

/* free memory allocated for basis set shell struct
 
   b - pointer to basis set shell array
   n - number of the shell */
struct basis_shell* basis_shell_free(struct basis_shell *b, unsigned n) {
  unsigned i;
  if (b) {
    for (i=0; i<n; i++) {
      b[i].exp = vec_ffree(b[i].exp);
      b[i].cf1 = vec_ffree(b[i].cf1);
      b[i].cf2 = vec_ffree(b[i].cf2);
      }
    free(b);
    }
  return(NULL);
  }

/* free memory allocated for basis set center struct
 
   b - pointer to basis set center array
   n - number of the center */
struct basis_center* basis_center_free(struct basis_center *b, unsigned n) {
  unsigned i;
  if (b) {
    for (i=0; i<n; i++)
      b[i].shell = basis_shell_free(b[i].shell,b[i].n_shells);
    free(b);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
