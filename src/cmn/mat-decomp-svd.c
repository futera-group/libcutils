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
#include "cmn/matrix.h"
#include "cmn/message.h"
#include "cmn/vector.h"

/* maximum number of iterations in single value decomposition */
#define MATH_SVD_MAX_ITER 10

/* -------------------------------------------------------------------------- */

/* single value decomposition - initial reduction to bidiagonal form
 
   u  - orthogonal n1*n2 matrix: u^t*u = i
   w  - diagonal elements of n2*n2 matrix
   n1 - number of rows in the a matrix
   n2 - number of columns in the a matrix
   a  - norm of the matrix
   f  - scaling factor
   g  - scaling factor
   t  - auxiliary working array */
void mat_fdecomp_svd_hauseholder(double **u, double *w, unsigned n1,
  unsigned n2, double *a, double *f, double *g, double *t) {
  double h,s,scale;
  unsigned i,j,k,l;
  (*g) = scale = (*a) = 0.0;
  for (i=0; i<n2; i++) {
    l = i+1;
    t[i] = scale*(*g);
    (*g) = s = scale = 0.0;
    if (i<n1) {
      for (k=i; k<n1; k++)
        scale += fabs(u[k][i]);
      if (scale) {
        for (k=i; k<n1; k++) {
          u[k][i] /= scale;
          s += (u[k][i]*u[k][i]);
          }
        (*f) = u[i][i];
        (*g) = ((*f)<0.0 ? sqrt(s) : -sqrt(s));
        h = ((*f)*(*g)-s);
        u[i][i] = ((*f)-(*g));
        for (j=l; j<n2; j++) {
          for (s=0.0,k=i; k<n1; k++)
            s += (u[k][i]*u[k][j]);
          (*f) = (s/h);
          for (k=i; k<n1; k++)
            u[k][j] += ((*f)*u[k][i]);
          }
        for (k=i; k<n1; k++)
          u[k][i] *= scale;
        }
      }
    w[i] = (scale*(*g));
    (*g) = s = scale = 0.0;
    if (i<n1 && (i+1)<n2) {
      for (k=l; k<n2; k++)
        scale += fabs(u[i][k]);
      if (scale) {
        for (k=l; k<n2; k++) { 
          u[i][k] /= scale;
          s += (u[i][k]*u[i][k]);
          }
        (*f) = u[i][l];
        (*g) = ((*f)<0.0 ? sqrt(s) : -sqrt(s));
        h = ((*f)*(*g)-s);
        u[i][l] = ((*f)-(*g));
        for (k=l; k<n2; k++)
          t[k] = (u[i][k]/h);
        for (j=l; j<n1; j++) {
          for (s=0.0,k=l; k<n2; k++)
            s += (u[j][k]*u[i][k]);
          for (k=l; k<n2; k++)
            u[j][k] += (s*t[k]);
          }
        for (k=l; k<n2; k++)
          u[i][k] *= scale;
        }
      }
    (*a) = math_fmax((*a),(fabs(w[i])+fabs(t[i])));
    }
  }

/* singular value decomposition of general matrix
 
   a  - the n1xn2 matrix that is to be decomposed: a = u*w*v^t
   u  - orthogonal n1*n2 matrix: u^t*u = i
   v  - orthogonal n2*n2 matrix: v^t*v = i
   w  - diagonal elements of n2*n2 matrix
   n1 - number of rows in the a matrix
   n2 - number of columns in the a matrix */
void mat_fdecomp_svd(double **a, double **u, double **v, double *w,
  unsigned n1, unsigned n2) {
  double c,f,g,h,s,x,y,z,anorm,*rv1;
  unsigned i,j,jj,k,l,nm,its;
  short flag;
  /* initialization */
  rv1 = vec_falloc(n2);
  mat_fcopy(u,a,n1,n2);
  /* householder reduction to bidiagonal form */
  mat_fdecomp_svd_hauseholder(u,w,n1,n2,&anorm,&f,&g,rv1);
  /* accumulation of right-hand transformation */
  for (i=n2; i>=1; i--) {
    if (i<n2) {
      if (g) {
        for (j=l; j<n2; j++)
          v[j][i-1] = ((u[i-1][j]/u[i-1][l])/g);
        for (j=l; j<n2; j++) {
          for (s=0.0,k=l; k<n2; k++)
            s += (u[i-1][k]*v[k][j]);
          for (k=l; k<n2; k++)
            v[k][j] += (s*v[k][i-1]);
          }
        }
      for (j=l; j<n2; j++)
        v[i-1][j] = v[j][i-1]=0.0;
      }
    v[i-1][i-1] = 1.0;
    g = rv1[i-1];
    l = (i-1);
    }
  /* accumulation of left-hand transformation */
  for (i=math_umin(n1,n2); i>=1; i--) {
    l = i;
    g = w[i-1];
    for (j=l; j<n2; j++)
      u[i-1][j] = 0.0;
    if (g) {
      g = (1.0/g);
      for (j=l; j<n2; j++) {
        for (s=0.0,k=l; k<n1; k++)
          s += (u[k][i-1]*u[k][j]);
        f = ((s/u[i-1][i-1])*g);
        for (k=(i-1); k<n1; k++)
          u[k][j] += (f*u[k][i-1]);
        }
      for (j=(i-1); j<n1; j++)
        u[j][i-1] *= g;
      }
    else {
      for (j=(i-1); j<n1; j++)
        u[j][i-1] = 0.0;
      }
    u[i-1][i-1]++;
    }
  /* diagonalization of the bidiagonal form */
  for (k=n2; k>=1; k--) {
    for (its=0; its<MATH_SVD_MAX_ITER; its++) {
      flag = 1;
      /* test for splitting */
      for (l=k; l>=1; l--) {
        nm = l-1;
        if (((double)(fabs(rv1[l-1])+anorm))==anorm) {
          flag = 0;
          break;
          }
        if (!nm || (((double)(fabs(w[nm-1])+anorm))==anorm))
          break;
        }
      /* cancellation of rv1[l] if l>0 */
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i=(l-1); i<k; i++) {
          f = (s*rv1[i]);
          rv1[i] = (c*rv1[i]);
          if (((double)(fabs(f)+anorm))==anorm)
            break;
          g = w[i];
          h = math_fpyth(f,g);
          w[i] = h;
          h = (1.0/h);
          c = (g*h);
          s = (-f*h);
          for (j=0; j<n1; j++) {
            y = u[j][nm-1];
            z = u[j][i];
            u[j][nm-1] = (y*c+z*s);
            u[j][i] = (z*c-y*s);
            }
          }
        }
      z = w[k-1];
      /* convergence */
      if (l==k) {
        if (z<0.0) {
          w[k-1] = (-z);
          for (j=0; j<n2; j++)
            v[j][k-1] = (-v[j][k-1]);
          }
        break;
        }
      if ((its+1)==MATH_SVD_MAX_ITER)
        msg_error_f("maximum number of iterations (%d) in svd exceeded",1, 
          MATH_SVD_MAX_ITER);
      /* shift from bottom 2x2 minor */
      x = w[l-1];
      nm = (k-1);
      y = w[nm-1];
      g = rv1[nm-1];
      h = rv1[k-1];
      f = ((((y-z)*(y+z))+((g-h)*(g+h)))/(2.0*h*y));
      g = math_fpyth(f,1.0);
      f = ((((x-z)*(x+z))+(h*((y/(f+(f<0.0 ? -g : g)))-h)))/x);
      /* next QR transformation */
      c = s = 1.0;
      for (j=(l-1); j<nm; j++) {
        i = (j+1);
        g = rv1[i];
        y = w[i];
        h = (s*g);
        g = (c*g);
        z = math_fpyth(f,h);
        rv1[j] = z;
        c = (f/z);
        s = (h/z);
        f = ((x*c)+(g*s));
        g = ((g*c)-(x*s));
        h = (y*s);
        y *= c;
        for (jj=0; jj<n2; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = ((x*c)+(z*s));
          v[jj][i] = ((z*c)-(x*s));
          }
        z = math_fpyth(f,h);
        w[j] = z;
        if (z) {
          z = (1.0/z);
          c = (f*z);
          s = (h*z);
          }
        f = ((c*g)+(s*y));
        x = ((c*y)-(s*g));
        for (jj=0; jj<n1; jj++) {
          y = u[jj][j];
          z = u[jj][i];
          u[jj][j] = ((y*c)+(z*s));
          u[jj][i] = ((z*c)-(y*s));
          }
        }
      rv1[l-1] = 0.0;
      rv1[k-1] = f;
      w[k-1] = x;
      }
    }
  /* clean memory */
  vec_ffree(rv1);
  }

/* -------------------------------------------------------------------------- */
