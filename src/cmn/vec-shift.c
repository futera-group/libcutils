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

/* -------------------------------------------------------------------------- */

/* shift vector elements right (unsigned int)

   v - unsigned integer vector
   n - number of elements in the vector */
void vec_urshift(unsigned *v, unsigned n) {
  unsigned i;
  for (i=n-1; i>0; i--)
    v[i] = v[i-1];
  v[0] = 0.0;
  }

/* insert value and shift vector elements right (unsigned int)

   v - unsigned integer vector
   n - number of elements in the vector
   r - element ID for insertion
   x - value for insertion */
void vec_urshift_id(unsigned *v, unsigned n, unsigned r, unsigned x) {
  unsigned i;
  /* shift vector */
  for (i=n-1; i>r; i--)
    v[i] = v[i-1];
  /* insert value */
  v[r] = x;
  }

/* -------------------------------------------------------------------------- */

/* shift vector elements left (double precision real)

   v - double precision real vector
   n - number of elements in the vector */
void vec_flshift(double *v, unsigned n) {
  unsigned i;
  for (i=0; i<n-1; i++)
    v[i] = v[i+1];
  v[n-1] = 0.0;
  }

/* shift vector elements right (double precision real)

   v - double precision real vector
   n - number of elements in the vector */
void vec_frshift(double *v, unsigned n) {
  unsigned i;
  for (i=n-1; i>0; i--)
    v[i] = v[i-1];
  v[0] = 0.0;
  }

/* insert value and shift vector elements right (double)

   v - double precision real vector
   n - number of elements in the vector
   r - element ID for insertion
   x - value for insertion */
void vec_frshift_id(double *v, unsigned n, unsigned r, double x) {
  unsigned i;
  /* shift vector */
  for (i=n-1; i>r; i--)
    v[i] = v[i-1];
  /* insert value */
  v[r] = x;
  }

/* -------------------------------------------------------------------------- */
