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

#include <cmn/quaternion.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* Calculate bond vector for selected atom
 
   r     - the bond vector (atom-anchor)
   c2    - first connection coordinates (bond)
   c3    - second connection coordinates (angle)
   c4    - third connection coordinates (dihedral)
   bond  - bond length
   angle - valence angle
   dihed - dihedral angle */
void mol_zmt_atom_vec(double *r, double *c2, double *c3, double *c4,
  double bond, double angle, double dihed) {
  double ab[3],bc[3],s[3],*qt;
  /* initialization */
  qt = quat_alloc();
  vec_fset(r,0.0,3);
  /* bond angles */
  vec_fsub(ab,c4,c3,3);
  vec_fsub(bc,c3,c2,3);
  vec_funit(ab,3);
  vec_funit(bc,3);
  vec_fprod_vector(ab,bc,s);
  vec_funit(s,3);
  /* angle rotation */
  vec_fadd_scaled(r,bc,-1.0,3);
  quat_from_axis(qt,s,180.0-angle);
  quat_vec_rot(qt,r,1);
  /* dihedral rotation */
  quat_from_axis(qt,bc,-dihed);
  quat_vec_rot(qt,r,1);
  /* bond length */
  vec_fscale(r,bond,3);
  /* clean memory */
  quat_free(qt);
  }

/* -------------------------------------------------------------------------- */
