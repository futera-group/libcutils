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
#include "cmn/message.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* integrate function by trapezoidal rule
 
   f - array with functional values (regular discretization)
   a - lower boundary of the integration
   w - bin width of the discretization
   n - number of bins */
double math_integ_trapz(double *f, double a, double w, unsigned n) {
  double s = 0.0;
  unsigned i;
  /* no points */
  if (n<1)
    s = 0.0;
  /* one point only */
  else if (n==1)
    s = (w*f[0]);
  /* two points only */
  else if (n==2)
    s = (0.5*w*(f[0]+f[1]));
  /* full array */
  else {
    s = ((f[0]+f[n-1])/2.0);
    for (i=1; i<(n-1); i++)
      s += f[i];
    s *= w;
    }
  return(s);
  }

/* integrate function by trapezoidal rule - 1/N3 accuracy
 
   f - array with functional values (regular discretization)
   a - lower boundary of the integration
   w - bin width of the discretization
   n - number of bins */
double math_integ_trapz3n(double *f, double a, double w, unsigned n) {
  double s=0.0;
  unsigned i;
  /* less than 5 points */
  if (n<5)
    s = math_integ_trapz(f,a,w,n);
  /* full array */
  else {
    s = ((5.0*(f[0]+f[n-1])+13.0*(f[1]+f[n-2]))/12.0);
    for (i=2; i<(n-2); i++)
      s += f[i];
    s *= w;
    }
  return(s);
  }

/* integrate function by trapezoidal rule - 1/N4 accuracy
 
   f - array with functional values (regular discretization)
   a - lower boundary of the integration
   w - bin width of the discretization
   n - number of bins */
double math_integ_trapz4n(double *f, double a, double w, unsigned n) {
  double s=0.0;
  unsigned i;
  /* less than 7 points */
  if (n<7)
    s = math_integ_trapz3n(f,a,w,n);
  /* full array */
  else {
    s = ((9.0*(f[0]+f[n-1])+28.0*(f[1]+f[n-2])+23.0*(f[2]+f[n-3]))/24.0);
    for (i=3; i<(n-3); i++)
      s += f[i];
    s *= w;
    }
  return(s);
  }

/* integrate function by Simpson's rule
 
   f - array with functional values (regular discretization)
   a - lower boundary of the integration
   w - bin width of the discretization
   n - number of bins */
double math_integ_simpson(double *f, double a, double w, unsigned n) {
  double s=0.0;
  unsigned i;
  /* less than 3 points */
  if (n<3)
    s = math_integ_trapz(f,a,w,n);
  /* full array */
  else {
    s = (f[0]+f[n-1]);
    for (i=1; i<(n-1); i+=2)
      s += (4.0*f[i]);
    for (i=2; i<(n-1); i+=2)
      s += (2.0*f[i]);
    s *= (w/3.0);
    }
  return(s);
  }

/* calculate spherical integral from given function 
   (longitude/latitude grid with centroid rule)
 
   r - radius of the sphere
   h - maximal length of triangle on the sphere
   f - integrand function
   d - data for the integrand function */
double math_integ_sphere_llc (double r, double h,
  double f(double*, void*), void *d) {
  double result,sphere_area,sector_area,area,r2;
  double theta,theta1,theta2,phi,phi1,phi2;
  double x[3],x11[3],x12[3],x21[3],x22[3],val;
  unsigned ip,it,phi_num,theta_num;
  if (r==0)
    return(0.0);
  if (r<=0.0)
    msg_error("positive sphere radius expected"
      " in math_integ_sphere_llc function",1);
  if (h<=0.0)
    msg_error("positive triangle side limit expected"
      " in math_integ_sphere_llc function",1);
  /* intialization */
  r2 = r*r;
  result = 0.0;
  /* choose phi and theta counts that make short sides */
  phi_num = (unsigned)(M_PI*r/h);
  if ((phi_num*h)<M_PI)
    phi_num++;
  theta_num = (unsigned)(2.0*M_PI*r/h);
  if ((theta_num*h)<M_PI)
    theta_num++;
  /* only one theta (and hence, only one phi) */
  if (theta_num==1) {
    sphere_area = 4.0*M_PI*r2;
    theta = 0.0;
    phi = M_PI/2.0;
    math_rpt_xyz(x,r,phi,theta);
    val = f(x,d);
    result = sphere_area*val;
    }
  /* more theta points for one phi angle */
  else if (phi_num==1) {
    sphere_area = 4.0*M_PI*r2;
    sector_area = (sphere_area/theta_num);
    for (it=0; it<theta_num; it++) {
      theta = (2*it*M_PI)/theta_num;
      phi = M_PI/2.0;
      math_rpt_xyz(x,r,phi,theta);
      val = f(x,d);
      result += (sector_area*val);
      }
    }
  /* general case */
  else {
    /* top row triangle below north pole */
    phi1 = 0.0;
    phi2 = M_PI/phi_num;
    for (it=0; it<theta_num; it++) {
      theta1 = (2.0*it*M_PI/theta_num);
      theta2 = (2.0*(it+1)*M_PI /theta_num);
      math_rpt_xyz(x11,1.0,phi1,theta1);
      math_rpt_xyz(x12,1.0,phi2,theta1);
      math_rpt_xyz(x22,1.0,phi2,theta2);
      area = r2*math_rpt_area_3p(x11,x12,x22);
      math_rpt_area_3p_center(x,x11,x12,x22);
      vec_fscale(x,r,3);
      val = f(x,d);
      result += (area*val);
      }
    /* intermediate rows with squares divided to two triangles */
    for (ip=1; ip<(phi_num-1); ip++) {
      phi1 = (M_PI*ip/phi_num);
      phi2 = (M_PI*(ip+1)/phi_num);
      for (it=0; it<theta_num; it++) {
        theta1 = (2.0*it*M_PI/theta_num);
        theta2 = (2.0*(it+1)*M_PI/theta_num);
        math_rpt_xyz(x11,1.0,phi1,theta1);
        math_rpt_xyz(x21,1.0,phi1,theta2);
        math_rpt_xyz(x12,1.0,phi2,theta1);
        math_rpt_xyz(x22,1.0,phi2,theta2);
        /* first triangle */
        area = r2*math_rpt_area_3p(x11,x12,x22);
        math_rpt_area_3p_center(x,x11,x12,x22);
        vec_fscale(x,r,3);
        val = f(x,d);
        result += (area*val);
        /* second triangle */
        area = r2*math_rpt_area_3p(x22,x21,x11);
        math_rpt_area_3p_center(x,x22,x21,x11);
        vec_fscale(x,r,3);
        val = f(x,d);
        result += (area*val);
        }
      }
    /* bottom row triangle above south pole */
    phi1 = ((phi_num-1)*M_PI/phi_num);
    phi2 = M_PI;
    for (it=0; it<theta_num; it++) {
      theta1 = (2.0*it*M_PI/theta_num);
      theta2 = (2.0*(it+1)*M_PI/theta_num);
      math_rpt_xyz(x11,1.0,phi1,theta1);
      math_rpt_xyz(x21,1.0,phi1,theta2);
      math_rpt_xyz(x22,1.0,phi2,theta2);
      area = r2*math_rpt_area_3p(x11,x22,x21);
      math_rpt_area_3p_center(x,x11,x22,x21);
      vec_fscale(x,r,3);
      val = f(x,d);
      result += (area*val);
      }
    }
  return(result);
  }

/* -------------------------------------------------------------------------- */
