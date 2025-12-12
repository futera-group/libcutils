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

#ifndef ZF_LIB_CMN_QUATERNION_H
#define ZF_LIB_CMN_QUATERNION_H

/* -------------------------------------------------------------------------- */

/* allocate memory for quaternion */
double *quat_alloc(void);
/* free memory allocated for quaternion */
void *quat_free(double*);

/* normalize quaternion to unit size */
void quat_norm(double*);

/* conjugate quaternion */
void quat_conj(double*);
/* new conjugated quaternion */
void quat_conj_c(double*, double*);

/* multiply two quaternions q = a*b */
void quat_mult(double*, double*, double*);

/* create quaternion from axis and rotation angle */
void quat_from_axis(double*, double*, double);
/* create quaternion from vector */
void quat_from_vec(double*, double*);
/* create quaternion from Euler angles */
void quat_from_euler(double*, double, double, double);

/* convert quaternion to axis and rotation angle */
void quat_to_angle(double*, double*, double*);
/* convert quaternion to vector */
void quat_to_vec(double*, double*);
/* convert quaternion to matrix */
void quat_to_mat(double*, double**);
/* convert quaternion to 3x3 matrix */
void quat_to_mat3(double*, double**);
/* convert quaternion to matrix (vector of columns) */
void quat_to_mat4c(double*, double*);
/* convert quaternion to matrix (vector of rows) */
void quat_to_mat4r(double*, double*);

/* use quaternion to rotate vector */
void quat_vec_rot(double*, double*, short);

/* -------------------------------------------------------------------------- */

#endif
