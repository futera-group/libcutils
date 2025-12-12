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
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* allocate 2D unsigned integer array with variable lengths

   n - number of rows in the array
   m - lengths of rows in the array */
unsigned **array_u2Valloc(unsigned n, unsigned *m) {
  unsigned i,**x = NULL;
  if (n) {
    x = (unsigned**)malloc(n*sizeof(unsigned*));
    if (!x)
      msg_error("cannot allocate new 2D variable length array",1);
    for (i=0; i<n; i++) {
      x[i] = NULL;
      x[i] = (unsigned*)malloc(m[i]*sizeof(unsigned));
      if (!x[i])
         msg_error("cannot allocate new 2D variable lenght array",1);
      }
    }
  return(x);
  }

/* allocated 3D unsigned integer array with fixed dimensions
 
   n1,n2,n3 - dimensions of the array */
unsigned ***array_u3alloc(unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k,***id = NULL;
  /* memory allocation */
  id = (unsigned***)malloc(n1*sizeof(unsigned**));
  if (!id)
    msg_error("cannot allocate memory for 3D u-int array",1);
  for (i=0; i<n1; i++) {
    id[i] = NULL;
    id[i] = (unsigned**)malloc(n2*sizeof(unsigned*));
    if (!(id[i]))
      msg_error("cannot allocate memory for 3D u-int array",1);
    for (j=0; j<n2; j++) {
      id[i][j] = NULL;
      id[i][j] = (unsigned*)malloc(n3*sizeof(unsigned));
      if (!(id[i][j]))
        msg_error("cannot allocate memory for 3D u-int array",1);
      }
    }
  /* initialization */
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        id[i][j][k] = 0;
  return(id);
  }

/* allocated 3D double-precision real array with fixed dimensions
 
   n1,n2,n3 - dimensions of the array */
double ***array_f3alloc(unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k;
  double ***id = NULL;
  /* memory allocation */
  id = (double***)malloc(n1*sizeof(double**));
  if (!id)
    msg_error("cannot allocate memory for 3D double array",1);
  for (i=0; i<n1; i++) {
    id[i] = NULL;
    id[i] = (double**)malloc(n2*sizeof(double*));
    if (!(id[i]))
      msg_error("cannot allocate memory for 3D double array",1);
    for (j=0; j<n2; j++) {
      id[i][j] = NULL;
      id[i][j] = (double*)malloc(n3*sizeof(double));
      if (!(id[i][j]))
        msg_error("cannot allocate memory for 3D double array",1);
      }
    }
  /* initialization */
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++)
        id[i][j][k] = 0.0;
  return(id);
  }

/* -------------------------------------------------------------------------- */

/* free memory of 2D unsigned integer array with variable lengths

   n - number of rows in the array */
void* array_u2Vfree(unsigned **a, unsigned n) {
  unsigned i;
  if (a!=NULL) {
    for (i=0; i<n; i++)
      free(a[i]);
    free(a);
    }
  return(NULL);
  }

/* free memory reserved for the 3D u-int array with fixed dimensions
 
   id - pointer to the array
   n1 - first array dimension
   n2 - second array dimension */
void array_u3free(unsigned ***id, unsigned n1, unsigned n2) {
  unsigned i,j;
  if (id) {
    for (i=0; i<n1; i++) {
      if (id[i]) {
        for (j=0; j<n2; j++)
          if (id[i][j])
            free(id[i][j]);
        }
      free(id[i]);
      }
    free(id);
    }
  }

/* free memory reserved for the 3D double array with fixed dimensions
 
   id - pointer to the array
   n1 - first array dimension
   n2 - second array dimension */
void array_f3free(double ***id, unsigned n1, unsigned n2) {
  unsigned i,j;
  if (id) {
    for (i=0; i<n1; i++) {
      if (id[i]) {
        for (j=0; j<n2; j++)
          if (id[i][j])
            free(id[i][j]);
        }
      free(id[i]);
      }
    free(id);
    }
  }

/* -------------------------------------------------------------------------- */

/* create copy of 2D unsigned integer array with variable lengths

   a - pointer to original array
   n - number of rows in the array
   m - lengths of rows in the array */
unsigned **array_u2Vcopy_new(unsigned **a, unsigned n, unsigned *m) {
  unsigned i,j,**x = NULL;
  x = array_u2Valloc(n,m);
  for (i=0; i<n; i++)
    for (j=0; j<m[i]; j++)
      x[i][j] = a[i][j];
  return(x);
  }

/* -------------------------------------------------------------------------- */

/* search given value in 2D unsigned integer array with variable lengths

   a - pointer to original array
   v - the value which is searched
   n - number of rows in the array
   m - lengths of rows in the array */
int array_u2Vfind(unsigned **a, unsigned v, unsigned n, unsigned *m) {
  unsigned i,j;
  for (i=0; i<n; i++)
    for (j=0; j<m[i]; j++)
      if (a[i][j]==v)
        return(i);
  return(-1);
  }

/* -------------------------------------------------------------------------- */

/* change number of rows in array with variable lengths
 
   a  - pointer to the array
   c  - number of columns in each row (after resizing)
   r1 - number of rows before resizing
   r2 - number of rows after resizing */
unsigned **array_u2Vresize_r(unsigned **a, unsigned *c, 
  unsigned r1, unsigned r2) {
  unsigned i,j,n,**b = NULL;
  if (r1==r2)
    b = a;
  else if (r2!=0) {
    b = array_u2Valloc(r2,c);
    n = (r1<r2 ? r1 : r2);
    for (i=0; i<n; i++)
      for (j=0; j<c[i]; j++)
        b[i][j] = a[i][j];
    array_u2Vfree(a,r1);
    }
  else
    array_u2Vfree(a,r1);
  return(b);
  }

/* -------------------------------------------------------------------------- */

/* delete specific row in array with variable lengths
 
   a - pointer to the array
   c - number of columns in each row (after deletion)
   r - number of rows before deletion
   d - ID of the row for deletion */
unsigned **array_u2Vdelete_r(unsigned **a, unsigned *c, 
  unsigned r, unsigned d) {
  unsigned i,j,**b = NULL;
  if (r>1) {
    b = array_u2Valloc(r-1,c);
    for (i=0; i<d; i++)
      for (j=0; j<c[i]; j++)
        b[i][j] = a[i][j];
    for (i=d; i<r-1; i++)
      for (j=0; j<c[i]; j++)
        b[i][j] = a[i+1][j];
    array_u2Vfree(a,r);
    }
  else
    array_u2Vfree(a,r);
  return(b);
  }

/* -------------------------------------------------------------------------- */
