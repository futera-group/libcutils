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

#include <stdlib.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for simulation cell data */
struct cell* cell_new(void) {
  struct cell* c = NULL;
  /* memory allocation */
  c = malloc(sizeof(struct cell));
  if (!c)
    msg_error("cannot allocate memory for simulation-cell data",1);
  /* initialization */
  c->type = 0;
  c->transf = 0;
  c->space_group_name = NULL;
  c->space_group_id = 0;
  c->origin = vec_falloc(3);
  vec_fset(c->origin,0.0,3);
  c->side = vec_falloc(3);
  vec_fset(c->side,0.0,3);
  c->angle = vec_falloc(3);
  vec_fset(c->angle,0.0,3);
  c->vector = mat_falloc(3,3);
  mat_fset(c->vector,0.0,3,3);
  c->volume = 0.0;
  c->t_rot_a = NULL;
  c->t_c2r_a = NULL;
  c->t_r2c_a = NULL;
  c->t_c2r_c = NULL;
  c->t_r2c_c = NULL;
  return(c);
  }

/* free memory allocated for simulation cell data
 
   c - the cell data */
struct cell* cell_free(struct cell *c) {
  if (c) {
    str_free(c->space_group_name);
    c->origin = vec_ffree(c->origin);
    c->side = vec_ffree(c->side);
    c->angle = vec_ffree(c->angle);
    c->vector = mat_ffree(c->vector,3);
    c->t_rot_a = mat_ffree(c->t_rot_a,3);
    c->t_c2r_a = mat_ffree(c->t_c2r_a,3);
    c->t_r2c_a = mat_ffree(c->t_r2c_a,3);
    c->t_c2r_c = mat_ffree(c->t_c2r_c,3);
    c->t_r2c_c = mat_ffree(c->t_r2c_c,3);
    free(c);
    c = NULL;
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */

/* copy a cell data structure

   c1 - the new cell
   c0 - the source cell */
void cell_copy(struct cell *c1, struct cell *c0) {
  /* cell type */
  c1->type = c0->type;
  c1->transf = c0->transf;
  c1->space_group_name = str_free(c1->space_group_name);
  c1->space_group_name = str_copy_new(c0->space_group_name);
  c1->space_group_id = c0->space_group_id;
  /* origin */
  if (!c1->origin)
    c1->origin = vec_falloc(3);
  vec_fcopy(c1->origin,c0->origin,3);
  /* side lengths */
  if (!c1->side)
    c1->side = vec_falloc(3);
  vec_fcopy(c1->side,c0->side,3);
  /* side angles */
  if (!c1->angle)
    c1->angle = vec_falloc(3);
  vec_fcopy(c1->angle,c0->angle,3);
  /* cell vectors */
  if (!c1->vector)
    c1->vector = mat_falloc(3,3);
  mat_fcopy(c1->vector,c0->vector,3,3);
  /* volume */
  c1->volume = c0->volume;
  /* rotation matrix */
  if (c0->t_rot_a && !c1->t_rot_a)
    c1->t_rot_a = mat_falloc(3,3);
  mat_fcopy(c1->t_rot_a,c0->t_rot_a,3,3);
  /* transformation matrices */
  if (c0->t_c2r_a && !c1->t_c2r_a)
    c1->t_c2r_a = mat_falloc(3,3);
  mat_fcopy(c1->t_c2r_a,c0->t_c2r_a,3,3);
  if (c0->t_r2c_a && !c1->t_r2c_a)
    c1->t_r2c_a = mat_falloc(3,3);
  mat_fcopy(c1->t_r2c_a,c0->t_r2c_a,3,3);
  if (c0->t_c2r_c && !c1->t_c2r_c)
    c1->t_c2r_c = mat_falloc(3,3);
  mat_fcopy(c1->t_c2r_c,c0->t_c2r_c,3,3);
  if (c0->t_r2c_c && !c1->t_r2c_c)
    c1->t_r2c_c = mat_falloc(3,3);
  mat_fcopy(c1->t_r2c_c,c0->t_r2c_c,3,3);
  }

/* create a new copy of the given cell data structure

   c - the cell data structure */
struct cell* cell_copy_new(struct cell *c) {
  struct cell *x;
  x = cell_new();
  cell_copy(x,c);
  return(x);
  }

/* -------------------------------------------------------------------------- */
