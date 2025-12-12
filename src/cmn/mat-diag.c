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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/vector.h>

#define DIAG_QL_IMPL_MAX_ITER 30

/* -------------------------------------------------------------------------- */

/* reduce symmetric matrix to tridiagonal form by Householder method

   m - input matrix which will replaced by orthogonal matrix for transformation
   d - vector for diagonal elements of final reduced matrix
   e - vector for off-diagonal elements of final reduced matrix
   n - dimension of the matrix */
void mat_diag_householder(double **m, double *d, double *e, unsigned n) {
  unsigned i,j,k;
  double scale,h,hh,g,f;
  /* n-2 subsequent transformations */
  for (i=n-1; i>0; i--) {
    h = scale = 0.0;
    if (i>1) {
      /* criterion for skipping transformation */
      for (k=0; k<i; k++)
        scale += fabs(m[i][k]);
      /* skip transformation */
      if (scale==0.0)
        e[i] = m[i][i-1];
      /* transform i-th matrix block */
      else {
        /* rescale i-th row of m and calculate sigma factor */
        for (k=0; k<i; k++) {
          m[i][k] /= scale;
          h += m[i][k]*m[i][k];
          }
        f = m[i][i-1];
        g = (f>=0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale*g;
        /* h=0.5*|u|^2 */
        h -= f*g;
        /* save vector u in i-th row of m */
        m[i][i-1]= f-g;
        f = 0.0;
        for (j=0; j<i; j++) {
          /* save u/h in i-th column of m */
          m[j][i] = m[i][j]/h;
          /* construct p=(m*u)/h */
          g = 0.0;
          for (k=0; k<=j; k++)
            g += m[j][k]*m[i][k];
          for (k=j+1; k<i; k++)
            g += m[k][j]*m[i][k];
          /* save element of p in e */
          e[j] = g/h;
          /* scalar product u.p */
          f += (e[j]*m[i][j]);
          }
        /* calculate K=(u.p)/(2h) */
        hh = f/(h+h);
        /* calculate q=p-Ku and save it in e */
        for (j=0; j<i; j++) {
          f = m[i][j];
          e[j] = g = e[j]-hh*f;
          for (k=0; k<=j; k++)
            m[j][k] -= (f*e[k]+g*m[i][k]);
          }
        }
      }
    else
      e[i] = m[i][i-1];
    d[i] = h;
    }
  /* create transformation matrix needed for eigenvector calculation */
  d[0] = 0.0;
  e[0] = 0.0;
  /* accumulate transformation matrices */
  for (i=0; i<n; i++) {
    /* skip first block */
    if (d[i]!=0.0) {
      for (j=0; j<i; j++) {
        g = 0.0;
        /* use u and u/h for p.q scalar product */
        for (k=0; k<i; k++)
          g += (m[i][k]*m[k][j]);
        for (k=0; k<i; k++)
          m[k][j] -= g*m[k][i];
        }
      }
    d[i] = m[i][i];
    m[i][i] = 1.0;
    /* set row and column to identity matrix */
    for (j=0; j<i; j++)
      m[j][i] = m[i][j] = 0.0;
    }
  }

/* diagonalize 3-diagonal matrix by QL algorithm with implicit shifts

   d - vector with diagonal elements of the matrix (eigenvalues on output)
   e - vector with off-diagonal elements of the matrix (e[0] arbitrary)
   z - matrix from Householder transformation or identity matrix (eigenvectors)
   n - dimension of the matrix */
void mat_diag_ql_implicit(double *d, double *e, double **z, unsigned n) {
  int i,j,k,m;
  unsigned iter;
  double s,r,p,g,f,dd,c,b,x;
  /* renumber elements of vector e */
  for (i=1; i<n; i++) 
    e[i-1] = e[i]; 
  e[n-1] = 0.0;
  /* iterative transformation */
  for (i=0; i<n; i++) {
    iter = 0;
    do {
      /* check if there is a single small subdiagonal element */
      for (m=i; m<n-1; m++) {
        dd = fabs(d[m])+fabs(d[m+1]);
        if ((fabs(e[m])+dd)==dd)
          break;
        }
      if (m!=i) {
        /* check number of iterations */
        if ((iter++)==DIAG_QL_IMPL_MAX_ITER)
          msg_error("maximum number of iteration in max_diag_ql_implicit" 
            " exceeded",1);
        /* calculate shift */
        g = (d[i+1]-d[i])/(2.0*e[i]);
        r = hypot(g,1.0);
        x = (g>=0.0 ? fabs(r) : -fabs(r));
        g = d[m]-d[i]+e[i]/(g+x);
        s = c = 1.0;
        p = 0.0;
        /* plane rotation of QL and Givens rotation */
        for (j=m-1; j>=i; j--) {
          f = s*e[j];
          b = c*e[j];
          e[j+1] = r = hypot(f,g);
          /* underflow */
          if (r==0.0) {
            d[j+1] -= p;
            e[m] = 0.0;
            break;
            }
          s = f/r;
          c = g/r;
          g = d[j+1]-p;
          r = (d[j]-g)*s+2.0*c*b;
          d[j+1] = g+(p = s*r);
          g = c*r-b;
          /* calculate eigenvectors */
          for (k=0; k<n; k++) {
            f = z[k][j+1];
            z[k][j+1] = s*z[k][j]+c*f;
            z[k][j] = c*z[k][j]-s*f;
            }
          }
        if (r==0.0 && j>1)
          continue;
        d[i] -= p;
        e[i] = g;
        e[m] = 0.0;
        }
      }
    while (m!=i);
    }
  }

/* sort eigenvalues into descending order by straiht insertion method

   d - vector of eigenvalues
   v - matrix of eigenvectors
   n - number of eigenvalues and eigenvectors */
void mat_diag_sort_dsc(double *d, double **v, unsigned n) {
  unsigned i,j,k;
  double p;
  for (i=0; i<n-1; i++) {
    p = d[k = i];
    for (j=i; j<n; j++)
      if (d[j]>p)
        p = d[k = j];
    if (k!=i) {
      d[k] = d[i];
      d[i] = p;
      for (j=0; j<n; j++) {
        p = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = p;
        }
      }
    }
  }

/* sort eigenvalues into ascending order by direct insertion method

   d - vector of eigenvalues
   v - matrix of eigenvectors
   n - number of eigenvalues and eigenvectors */
void mat_diag_sort_asc(double *d, double **v, unsigned n) {
  unsigned i,j,k;
  double p;
  for (i=0; i<n-1; i++) {
    p = d[k = i];
    for (j=i; j<n; j++)
      if (d[j]<=p)
        p = d[k = j];
    if (k!=i) {
      d[k] = d[i];
      d[i] = p;
      for (j=0; j<n; j++) {
        p = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = p;
        }
      }
    }
  }

/* diagonalize square matrix, calculate eigenvalues and eigenvectors

   m - matrix which should be diagonalized
   u - matrix for eigenvectors
   e - vector for eigenvalues
   n - dimenstion of the matrix
   s - specify sort of eigenvalues */
void mat_diag(double **m, double **u, double *e, unsigned n, short s) {
  double **a,*v;
  /* memory allocation */
  a = mat_fcopy_new(m,n,n);
  v = vec_falloc(n);
  /* initialization */
  mat_fset(u,0.0,n,n);
  vec_fset(e,0.0,n);
  vec_fset(v,0.0,n);
  /* transform matrix to tridiagonal form */
  mat_diag_householder(a,e,v,n);
  /* diagonalize tridiagonal matrix */
  mat_diag_ql_implicit(e,v,a,n);
  /* sort eigenvalues */
  switch (s) {
    case MAT_DIAG_SORT_ASC: mat_diag_sort_asc(e,a,n); break;
    case MAT_DIAG_SORT_DSC: mat_diag_sort_dsc(e,a,n); break;
    }
  /* copy eigenvectors */
  mat_fcopy(u,a,n,n);
  /* free allocated memory */
  mat_ffree(a,n);
  vec_ffree(v);
  }

/* -------------------------------------------------------------------------- */
