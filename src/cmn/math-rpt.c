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
#include "cmn/math.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* converts spherical coordinates (r,p,t) to cartesian (x,y,z)
 
   x - array for cartesian coordinates
   r - radial distance
   p - polar angle
   t - azimutal angle */
void math_rpt_xyz(double *x, double r, double t, double p) {
  x[0] = r*cos(p)*sin(t);
  x[1] = r*sin(p)*sin(t);
  x[2] = r*cos(t);
  }

/* calculate lenght of triangle sides on the unit sphere
 
   v1,v2,v3 - vertices of the triangle (input)
   as,bs,cs - side lengths (output) */
void math_rpt_area_3p_sides(double *v1, double *v2, double *v3,
  double *as, double *bs, double *cs) {
  (*as) = math_acos(vec_fprod_scalar(v2,v3,3));
  (*bs) = math_acos(vec_fprod_scalar(v3,v1,3));
  (*cs) = math_acos(vec_fprod_scalar(v1,v2,3));
  }

/* calculate angles between triangle sides on the unit sphere
 
   as,bs,cs - the triangle side lenghts (input)
   a,b,c    - values of the angles (output) */
void math_rpt_area_3p_angles(double as, double bs, double cs, 
  double *a, double *b, double *c) {
  double ss,tan_a2,tan_b2,tan_c2;
  ss = (as+bs+cs)/2.0;
  tan_a2 = sqrt((sin(ss-bs)*sin(ss-cs))/(sin(ss)*sin(ss-as)));
  (*a) = (2.0*atan(tan_a2));
  tan_b2 = sqrt((sin(ss-as)*sin(ss-cs))/(sin(ss)*sin(ss-bs)));
  (*b) = (2.0*atan(tan_b2));
  tan_c2 = sqrt((sin(ss-as)*sin(ss-bs))/(sin(ss)*sin(ss-cs)));
  (*c) = (2.0*atan(tan_c2));
  }

/* calculate are selected by 3 points on the unit sphere
 
   v1,v2,v3 - coordinates of the points */
double math_rpt_area_3p(double *v1, double *v2, double *v3) {
  double a,b,c,as,bs,cs;
  /* compute the lengths of the sides of the spherical triangle */
  math_rpt_area_3p_sides(v1,v2,v3,&as,&bs,&cs);
  /* get the spherical angles */
  math_rpt_area_3p_angles(as,bs,cs,&a,&b,&c);
  /* the triangle area */
  return(a+b+c-M_PI);
  }

/* calculate coordinates of spherical triangle centroid
 
   vs       - coordinates of the centroid (output)
   v1,v2,v3 - three triangle vertices (input) */
void math_rpt_area_3p_center(double *vs, double *v1, double *v2, double *v3) {
  unsigned i;
  for (i=0; i<3; i++)
    vs[i] = (v1[i]+v2[i]+v3[i])/3.0;
  vec_funit(vs,3);
  }

/* -------------------------------------------------------------------------- */
