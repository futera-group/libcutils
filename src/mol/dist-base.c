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

#include <math.h>

/* -------------------------------------------------------------------------- */

/* distance between two points

   p1,p2 - cartesian coordinates of the two points
   dim   - dimension of the space */
double dist_r(double *p1, double *p2, unsigned dim) {
  unsigned i;
  double dist=0.0;
  for (i=0; i<dim; i++)
    dist += pow(p1[i]-p2[i],2.0);
  return(sqrt(dist));
  }

/* squared distance between two points

   p1,p2 - cartesian coordinates of the two points
   dim   - dimension of the space */
double dist_r2(double *p1, double *p2, unsigned dim) {
  unsigned i;
  double dist=0.0;
  for (i=0; i<dim; i++)
    dist += pow(p1[i]-p2[i],2.0);
  return(dist);
  }

/* -------------------------------------------------------------------------- */
