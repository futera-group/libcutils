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

#ifndef ZF_LIB_MOL_CELL_H
#define ZF_LIB_MOL_CELL_H

/* -------------------------------------------------------------------------- */

/* box types */
#define CELL_TYPE_TRIC  1  /* triclinic */
#define CELL_TYPE_MONO  2  /* monoclinic */
#define CELL_TYPE_ORTH  3  /* orthorhombic */
#define CELL_TYPE_RHOM  4  /* rhombohedral */
#define CELL_TYPE_TETR  5  /* tetragonal */
#define CELL_TYPE_HEXA  6  /* hexagonal */
#define CELL_TYPE_CUBC  7  /* cubic */

#define CELL_CHECK_ACC  1.0e-5  /* accuracy */

/* -------------------------------------------------------------------------- */

/* atomic data */
struct cell {
  short type;                      /* cell type */
  short transf;                    /* coordinate-transformation indicator */
  char* space_group_name;          /* space group name */
  unsigned space_group_id;         /* space group number */
  double *origin;                  /* origin of the cell */
  double *side;                    /* side lengths (a,b,c) */
  double *angle;                   /* cell angles (alpha,beta,gamma) */
  double **vector;                 /* cell vectors (a1,a2,a3) */
  double volume;                   /* volume of the cell */
  double **t_rot_a;                /* transf. matrix: axes rotation */
  double **t_c2r_a;                /* transf. matrix: cell -> ref (axes) */
  double **t_r2c_a;                /* transf. matrix: ref -> cell (axes) */
  double **t_c2r_c;                /* transf. matrix: cell -> ref (coord) */
  double **t_r2c_c;                /* transf. matrix: ref -> cell (coord) */
  };

/* -------------------------------------------------------------------------- */

/* allocate memory for simulation cell data */
struct cell* cell_new(void);
/* free memory allocated for simulation cell data */
struct cell* cell_free(struct cell*);

/* copy a cell data structure */
void cell_copy(struct cell*, struct cell*);
/* create a new copy of the given cell data structure */
struct cell* cell_copy_new(struct cell*);

/* set cell parameters by given side lengths */
void cell_set_side(struct cell*, double, double, double);
/* set cell parameters by given array of side lengths */
void cell_set_side_v(struct cell*, double*);
/* set cell parameters by given side lengths and angles */
void cell_set_side_angle(struct cell*, double, double, double,
  double, double, double);
/* set cell parameters by given arrays of side lengths and angles */
void cell_set_side_angle_v(struct cell*, double*, double*);
/* set cell parameters by given side lengths and tilting factors */
void cell_set_side_tilt(struct cell*, double, double, double,
  double, double, double);
/* set cell parameters by given arrays of side lengths and tilting factors */
void cell_set_side_tilt_v(struct cell*, double*, double*);
/* set cell parameters by given cell vectors */
void cell_set_vec(struct cell*, double*, double*, double*);

/* calculate titlting factors of the cell planes */
void cell_tilt(struct cell*, double*, double*, double*,
  double*, double*, double*);
/* calculate arrays of titlting factors of the cell planes */
void cell_tilt_v(struct cell*, double*, double*);

/* return internal ID of cell type */
short cell_type_id(char*);
/* convert internal cell type ID to corresponding flag */
char* cell_type_name(short, char*);
/* deduce and set cell type from cell vectors */
void cell_type_set(struct cell*);

/* calculate axes-rotation matrix aligning A1 to X and placing A2 to XY plane */
void cell_rmat_set(struct cell*);

/* check if the coordinate transformation is needed and set the matrices */
void cell_tmat_set(struct cell*);

/* -------------------------------------------------------------------------- */

#endif
