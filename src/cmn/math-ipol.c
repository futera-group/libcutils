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
#include "cmn/message.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* initialization values */
#define INIT_MIN   9.9e+99  /* minimum search */
#define INIT_MAX  -9.9e+99  /* maximum search */

/* -------------------------------------------------------------------------- */

/* find interval to which belongs the interpolation point

   x  - the interpolation point
   v  - array of X values
   nv - number of the X values
   id - ID of interval */
short math_ipol_locate(double x, double *v, long unsigned n,
  long unsigned *id) {
  long unsigned i0,i1,im;
  short ascnd;
  /* out of range */
  if (x<v[0] || x>v[n-1])
    return(0);
  /* initialization */
  i0 = 0;
  i1 = n-1;
  ascnd = (v[n-1]>=v[0] ? 1 : 0);
  /* interval bisection */
  while ((i1-i0)>1) {
    /* midpoint */
    im = (i0+i1)>>1;
    /* update interval */
    if ((x>=v[im]) && ascnd)
      i0 = im;
    else
      i1 = im;
    }
  /* position */
  (*id) = i0;
  return(1);
  }

/* linear interpolation between two points

   x  - the X value in [x0,x1] interval
   x0 - lower boundary X coordinate
   y0 - lower boundary Y coordinate
   x1 - upper boundary X coordinate
   y1 - upper boundary Y coordinate */
double math_ipol_linear(double x, double x0, double y0, double x1, double y1) {
  double df;
  if (x0==x1)
    return(0.0);
  df = (y1-y0)/(x1-x0);
  return(y0+df*(x-x0));
  }

/* polynomial interpolation between two points

   x  - the X value 
   vx - array of X values
   vy - array of Y values
   nv - number of values */
double math_ipol_polynomial(double x, double *vx, double *vy,
  long unsigned nv) {
  double t,dx,*c,*d,hp,ho,dh,dv,w,v = 0.0;
  long unsigned i,j,id = 0;
  /* memory allocation */
  c = vec_falloc(nv);
  d = vec_falloc(nv);
  /* initialization */
  dx = fabs(x-vx[0]);
  for (i=0; i<nv; i++) {
    /* nearest input value */
    t = fabs(x-vx[i]);
    if (t<dx) {
      id = i;
      dx = t;
      }
    /* working arrays */
    c[i] = vy[i];
    d[i] = vy[i];
    }
  /* initial approximation */
  v = vy[id--];
  /* neville's recursive algorithm */
  for (j=1; j<nv; j++) {
    for (i=0; i<nv-j; i++) {
      ho = vx[i]-x;
      hp = vx[i+j]-x;
      dh = ho-hp;
      if (dh==0.0)
        msg_error("problem in polynomial interpolation"
          " (identical input values)",1);
      w = c[i+1]-d[i];
      dh = w/dh;
      d[i] = hp*dh;
      c[i] = ho*dh;
      }
    /* update interpolation */
    dv = 2.0*(id+1);
    v += (dv<(nv-j) ? c[id+1] : d[id--]);
    }
  /* clean memory */
  vec_ffree(c);
  vec_ffree(d);
  return(v);
  }

/* calculate second derivatives of cubic-spline at given points

   vx - array of X values
   vy - array of Y values
   nv - number of values
   y0 - lower boundary condition
   y1 - upper boundary condition */
double* math_ipol_spline_y2(double *vx, double *vy, long unsigned nv,
  double y0, double y1) {
  long unsigned i;
  double *y2,*u,sgn,p,qn,un;
  /* memory allocation */
  y2 = vec_falloc(nv);
  u = vec_falloc(nv-1);
  /* lower boudary conditions */
  if (y0>=INIT_MIN) {
    /* natural condition */
    y2[0] = 0.0;
    u[0] = 0.0;
    }
  else {
    /* specified first derivative */
    y2[0] = -0.5;
    u[0] = (3.0/(vx[1]-vx[0]))*((vy[1]-vy[0])/(vx[1]-vx[0])-y0);
    }
  /* decomposition loop of the tridiagonal algorithm */
  for (i=1; i<nv-1; i++) {
    sgn = (vx[i]-vx[i-1])/(vx[i+1]-vx[i-1]);
    p = sgn*y2[i-1] + 2.0;
    y2[i] = (sgn-1.0)/p;
    u[i] = (vy[i+1]-vy[i])/(vx[i+1]-vx[i]) - (vy[i]-vy[i-1])/(vx[i]-vx[i-1]);
    u[i] = (6.0*u[i]/(vx[i+1]-vx[i-1])-sgn*u[i-1])/p;
    }
  /* upper boudary conditions */
  if (y1>=INIT_MIN) {
    /* natural condition */
    qn = 0.0;
    un = 0.0;
    }
  else {
    /* specified first derivative */
    qn = 0.5;
    un = (3.0/(vx[nv-1]-vx[nv-2]))*(y1-(vy[nv-1]-vy[nv-2])/(vx[nv-1]-vx[nv-2]));
    }
  y2[nv-1] = (un-qn*u[nv-2])/(qn*y2[nv-2]+1.0);
  /* back-substitution of the tridiagonal algorithm */
  for (i=nv-2; (i+1)>=1; i--)
    y2[i] = y2[i]*y2[i+1] + u[i];
  /* clean memory */
  vec_ffree(u);
  return(y2);
  }

/* cubic-spline interpolation between two points

   x  - the X value 
   vx - array of X values
   vy - array of Y values
   nv - number of values
   id - interpolation interval ID */
double math_ipol_spline(double x, double *vx, double *vy,
  long unsigned nv, long unsigned id) {
  long unsigned k0,k1;
  double *y2,dx,a,b,v;
  /* initialization */
  k0 = id;
  k1 = id+1;
  dx = vx[k1]-vx[k0];
  if (dx==0.0)
    msg_error("invalid data for spline interpolation",1);
  a = (vx[k1]-x)/dx;
  b = (x-vx[k0])/dx;
  /* second derivatives */
  y2 = math_ipol_spline_y2(vx,vy,nv,INIT_MIN,INIT_MIN);
  /* interpolation */
  v = a*vy[k0] + b*vy[k1] + ((a*a*a-a)*y2[k0]+(b*b*b-b)*y2[k1])*(dx*dx)/6.0;
  /* clean memory */
  vec_ffree(y2);
  return(v);
  }

/* -------------------------------------------------------------------------- */
