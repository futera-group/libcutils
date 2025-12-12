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
#include "cmn/message.h"
#include "cmn/vector.h"

#define QUAT_ZERO 1.0E-10

/* -------------------------------------------------------------------------- */

/* allocate memory for quaternion */
double *quat_alloc(void) {
  double *q = NULL;
  q = calloc(4,sizeof(double));
  if (!q)
    msg_error("cannot allocate memory for quaternion",1);
  return(q);
  }

/* free memory allocated for quaternion

   q - the quaternion */
void* quat_free(double *q) {
  if (q)
    free(q);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */

/* normalize quaternion to unit size

   q - the quaternion */
void quat_norm(double *q) {
  double qs = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
  unsigned i;
  if (qs>QUAT_ZERO && fabs(qs-1.0)>QUAT_ZERO) {
    qs = sqrt(qs);
    for (i=0; i<4; i++)
      q[i] /= qs;
    }
  }

/* -------------------------------------------------------------------------- */

/* conjugate quaternion

   q - the quaternion */
void quat_conj(double *q) { 
  q[1] = -q[1];
  q[2] = -q[2];
  q[3] = -q[3];
  }

/* new conjugated quaternion
 
   qc - quaternion copy 
   q  - original quaternion */
void quat_conj_c(double *qc, double *q) { 
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
  }

/* -------------------------------------------------------------------------- */

/* multiply two quaternions q = a*b
 
   q   - resulting quaternion product
   a,b - multipliers */
void quat_mult(double *q, double *a, double *b) {
  q[0] = a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
  q[1] = a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2];
  q[2] = a[0]*b[2]+a[2]*b[0]+a[3]*b[1]-a[1]*b[3];
  q[3] = a[0]*b[3]+a[3]*b[0]+a[1]*b[2]-a[2]*b[1];
  }

/* -------------------------------------------------------------------------- */

/* create quaternion from vector

   q - resulting quaternion
   v - the vector */
void quat_from_vec(double *q, double *v) {
  q[0] = 0.0;
  q[1] = v[0];
  q[2] = v[1];
  q[3] = v[2];
  }

/* create quaternion from axis a and rotation angle t

   q - the quaternion
   a - the axis
   t - the vector */
void quat_from_axis(double *q, double *a, double t) {
  double rd=M_PI/180.0;
  double st=sin(0.5*t*rd);
  vec_funit(a,3);
  q[0] = cos(0.5*t*rd);
  q[1] = a[0]*st;
  q[2] = a[1]*st;
  q[3] = a[2]*st;
  }

/* create quaternion from Euler angles

   q     - the quaternion
   a,b,c - the Euler angles */
void quat_from_euler(double *q, double a, double b, double c) {
  double sa,sb,sc,ca,cb,cc;
  double r,rd = M_PI/360.0;
  r = a*rd; sa = sin(r); ca = cos(r);
  r = b*rd; sb = sin(r); cb = cos(r);
  r = c*rd; sc = sin(r); cc = cos(r);
  q[0] = cc*ca*cb+sc*sa*sb;
  q[1] = sc*ca*cb-cc*sa*sb;
  q[2] = cc*sa*cb+sc*ca*sb;
  q[3] = cc*ca*sb-sc*sa*cb;
  quat_norm(q);
  }

/* -------------------------------------------------------------------------- */

/* convert quaternion to rotation axis and angle

   q - the quaternion
   a - the axis
   t - the angle */
void quat_to_angle(double *q, double *a, double *t) {
  double sc=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  a[0] = q[1]/sc;
  a[1] = q[2]/sc;
  a[2] = q[3]/sc;
  (*t) = 2.0*acos(q[0]);
  }

/* convert quaternion to vector

   q - the quaternion
   v - resulting vector */
void quat_to_vec(double *q, double *v) {
  v[0] = q[1];
  v[1] = q[2];
  v[2] = q[3];
  }

/* convert quaternion to 4x4 matrix

   q - the quaternion
   m - resulting matrix */
void quat_to_mat(double *q, double **m) {
  double x2,y2,z2,xy,xz,yz,wx,wy,wz;
  x2 = q[1]*q[1];
  y2 = q[2]*q[2];
  z2 = q[3]*q[3];
  xy = q[1]*q[2];
  xz = q[1]*q[3];
  yz = q[2]*q[3];
  wx = q[0]*q[1];
  wy = q[0]*q[2];
  wz = q[0]*q[3];
  m[0][0] = 1.0-2.0*(y2+z2);
  m[0][1] = 2.0*(xy-wz);
  m[0][2] = 2.0*(xz+wy);
  m[0][3] = 0.0;
  m[1][0] = 2.0*(xy+wz);
  m[1][1] = 1.0-2.0*(x2+z2);
  m[1][2] = 2.0*(yz-wx);
  m[1][3] = 0.0;
  m[2][0] = 2.0*(xz-wy);
  m[2][1] = 2.0*(yz+wx);
  m[2][2] = 1.0-2.0*(x2+y2);
  m[2][3] = 0.0;
  m[3][0] = 0.0;
  m[3][1] = 0.0;
  m[3][2] = 0.0;
  m[3][3] = 1.0;
  }

/* convert quaternion to 3x3 matrix

   q - the quaternion
   m - resulting matrix */
void quat_to_mat3(double *q, double **m) {
  double x2,y2,z2,xy,xz,yz,wx,wy,wz;
  x2 = q[1]*q[1];
  y2 = q[2]*q[2];
  z2 = q[3]*q[3];
  xy = q[1]*q[2];
  xz = q[1]*q[3];
  yz = q[2]*q[3];
  wx = q[0]*q[1];
  wy = q[0]*q[2];
  wz = q[0]*q[3];
  m[0][0] = 1.0-2.0*(y2+z2);
  m[0][1] = 2.0*(xy-wz);
  m[0][2] = 2.0*(xz+wy);
  m[1][0] = 2.0*(xy+wz);
  m[1][1] = 1.0-2.0*(x2+z2);
  m[1][2] = 2.0*(yz-wx);
  m[2][0] = 2.0*(xz-wy);
  m[2][1] = 2.0*(yz+wx);
  m[2][2] = 1.0-2.0*(x2+y2);
  }

/* convert quaternion to matrix (vector of columns)

   q - the quaternion
   m - resulting matrix */
void quat_to_mat4c(double *q, double *m) {
  double x2,y2,z2,xy,xz,yz,wx,wy,wz;
  x2 = q[1]*q[1];
  y2 = q[2]*q[2];
  z2 = q[3]*q[3];
  xy = q[1]*q[2];
  xz = q[1]*q[3];
  yz = q[2]*q[3];
  wx = q[0]*q[1];
  wy = q[0]*q[2];
  wz = q[0]*q[3];
  m[ 0] = 1.0-2.0*(y2+z2);
  m[ 1] = 2.0*(xy+wz);
  m[ 2] = 2.0*(xz-wy);
  m[ 3] = 0.0;
  m[ 4] = 2.0*(xy-wz);
  m[ 5] = 1.0-2.0*(x2+z2);
  m[ 6] = 2.0*(yz+wx);
  m[ 7] = 0.0;
  m[ 8] = 2.0*(xz+wy);
  m[ 9] = 2.0*(yz-wx);
  m[10] = 1.0-2.0*(x2+y2);
  m[11] = 0.0;
  m[12] = 0.0;
  m[13] = 0.0;
  m[14] = 0.0;
  m[15] = 1.0;
  }

/* convert quaternion to matrix (vector of rows)

   q - the quaternion
   m - resulting matrix */
void quat_to_mat4r(double *q, double *m) {
  double x2,y2,z2,xy,xz,yz,wx,wy,wz;
  x2 = q[1]*q[1];
  y2 = q[2]*q[2];
  z2 = q[3]*q[3];
  xy = q[1]*q[2];
  xz = q[1]*q[3];
  yz = q[2]*q[3];
  wx = q[0]*q[1];
  wy = q[0]*q[2];
  wz = q[0]*q[3];
  m[ 0] = 1.0-2.0*(y2+z2);
  m[ 1] = 2.0*(xy-wz);
  m[ 2] = 2.0*(xz+wy);
  m[ 3] = 0.0;
  m[ 4] = 2.0*(xy+wz);
  m[ 5] = 1.0-2.0*(x2+z2);
  m[ 6] = 2.0*(yz-wx);
  m[ 7] = 0.0;
  m[ 8] = 2.0*(xz-wy);
  m[ 9] = 2.0*(yz+wx);
  m[10] = 1.0-2.0*(x2+y2);
  m[11] = 0.0;
  m[12] = 0.0;
  m[13] = 0.0;
  m[14] = 0.0;
  m[15] = 1.0;
  }

/* -------------------------------------------------------------------------- */

/* use quaternion to rotate vector

   q - the quaternion
   v - the vector
   n - normalize the vector */
void quat_vec_rot(double *q, double *v, short n) {
  double vq[4],cq[4],rq[4],sq[4];
  if (n)
    vec_funit(v,3);
  quat_from_vec(vq,v);
  quat_conj_c(cq,q);
  quat_mult(rq,vq,cq);
  quat_mult(sq,q,rq);
  quat_to_vec(sq,v);
  }

/* -------------------------------------------------------------------------- */
