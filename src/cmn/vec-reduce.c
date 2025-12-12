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

#include <string.h>
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* reduce vector according to given mask (integer version)

   v    - data vector which will be reduced
   vn   - number of elements in the data vector
   m    - mask vector which IDs of element which will be selected
   mn   - number of elements in the mask vector */
int *vec_ireduce_mask(int *v, unsigned vn, unsigned *m, unsigned mn) {
  unsigned i,j=0;
  int *nv = vec_ialloc(mn);
  for (i=0; i<vn; i++)
    if (vec_ufind(m,i,mn)>=0)
      nv[j++] = v[i];
  vec_ifree(v);
  return(nv);
  }

/* reduce vector to given range (integer version)

   v  - data vector which will be reduced
   n  - number of elements in the data vector
   r0 - ID of the first element in the range
   r1 - ID of the last element in the range */
int *vec_ireduce_range(int *v, unsigned n, unsigned r0, unsigned r1) {
  unsigned i,nn = r1-r0+1;
  int *w = vec_ialloc(nn);
  for (i=0; i<nn; i++)
    w[i] = v[r0+i];
  vec_ifree(v);
  return(w);
  }

/* -------------------------------------------------------------------------- */

/* delete multiple elements and leave only one value of it (u-int)
 
   v - vector of unsigned integers
   n - number of elements in the vector (will be updated) */
unsigned *vec_ureduce_unique(unsigned *v, unsigned *n) {
  unsigned i,j,id=0,m=(*n),*vn;
  short null;
  /* empty vector or one element only */
  if (!v || m<2)
    return(v);
  /* find out if there is zero */
  null = (vec_ufind(v,0,*n)<0 ? 0 : 1);
  /* set non-unique elements to zero */
  for (i=0; i<(*n)-1; i++) {
    if (v[i]>0) {
      for (j=i+1; j<(*n); j++)
        if (v[i]==v[j]) {
          v[j]=0;
          m--;
          }
      }
    }
  /* no non-unique elements found */
  if ((*n)==m)
    return(v);
  /* copy remaining values to new vector */
  vn = vec_ualloc(m);
  for (i=0; i<(*n); i++)
    if (v[i] || null) {
      if (!v[i])
        null = 0;
      vn[id] = v[i];
      id++;
      }
  (*n) = m;
  /* free original vector */
  vec_ufree(v);
  return(vn);
  }

/* -------------------------------------------------------------------------- */

/* reduce vector according to given mask (double precision real version)

   v    - data vector which will be reduced
   vn   - number of elements in the data vector
   m    - mask vector which IDs of element which will be selected
   mn   - number of elements in the mask vector
   smax - lenght of single string */
double *vec_freduce_mask(double *v, unsigned vn, unsigned *m, unsigned mn) {
  unsigned i,j = 0;
  double *nv = vec_falloc(mn);
  for (i=0; i<vn; i++)
    if (vec_ufind(m,i,mn)>=0)
      nv[j++] = v[i];
  vec_ffree(v);
  return(nv);
  }

/* reduce vector to given range (double precision real version)

   v  - data vector which will be reduced
   n  - number of elements in the data vector
   r0 - ID of the first element in the range
   r1 - ID of the last element in the range */
double *vec_freduce_range(double *v, unsigned n, unsigned r0, unsigned r1) {
  unsigned i,nn = r1-r0+1;
  double *w = vec_falloc(nn);
  for (i=0; i<nn; i++)
    w[i] = v[r0+i];
  vec_ffree(v);
  return(w);
  }

/* -------------------------------------------------------------------------- */

/* reduce vector according to given mask (string version)

   v    - data vector which will be reduced
   vn   - number of elements in the data vector
   m    - mask vector which IDs of element which will be selected
   mn   - number of elements in the mask vector
   smax - lenght of single string */
char **vec_sreduce_mask(char **v, unsigned vn, unsigned *m, unsigned mn,
  unsigned smax) {
  unsigned i,j = 0;
  char **nv = vec_salloc(mn,smax);
  for (i=0; i<vn; i++)
    if (vec_ufind(m,i,mn)>=0)
      strcpy(nv[j++],v[i]);
  vec_sfree(v,vn);
  return(nv);
  }

/* reduce vector to given range (string version)

   v  - data vector which will be reduced
   n  - number of elements in the data vector
   r0 - ID of the first element in the range
   r1 - ID of the last element in the range
   m  - lenght of strings in the array */
char **vec_sreduce_range(char **v, unsigned n, unsigned r0, unsigned r1,
  unsigned m) {
  unsigned i,nn = r1-r0+1;
  char **w = vec_salloc(nn,m);
  for (i=0; i<nn; i++) 
    strncpy(w[i],v[r0+i],m);
  vec_sfree(v,n);
  return(w);
  }

/* -------------------------------------------------------------------------- */
