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
#include <stdlib.h>
#include <time.h>

/* -------------------------------------------------------------------------- */

/* set initial seed for random number generator

   s - the seed (if less than zero, actual time is used as seed) */
void rnd_init(int s) {
  if (s<0)
    srand(time(NULL));
  else
    srand(s);
  }

/* return random number from uniform distribution */
double rnd_uniform(void) {
  return(((double)rand())/RAND_MAX);
  }

/* return random number from Gaussian distribution */
double rnd_gauss(void) {
  double x1,x2,w;
  do {
    x1 = 2.0*((double)rand())/RAND_MAX-1.0;
    x2 = 2.0*((double)rand())/RAND_MAX-1.0;
    w = x1*x1+x2*x2;
    }
  while (w>=1.0);
  w = sqrt((-2.0*log(w))/w);
  if (((double)rand())/RAND_MAX<0.5)
    return(x1*w);
  return(x2*w);
  }

/* -------------------------------------------------------------------------- */
