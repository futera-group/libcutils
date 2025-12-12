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
#include <stdio.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of cell type
 
   flag - the type flag */
short cell_type_id(char *flag) {
  if (str_compare(flag,"triclinic"))
    return(CELL_TYPE_TRIC);
  else if (str_compare(flag,"monoclinic"))
    return(CELL_TYPE_MONO);
  else if (str_compare(flag,"orthorhombic"))
    return(CELL_TYPE_ORTH);
  else if (str_compare(flag,"rhombohedral"))
    return(CELL_TYPE_RHOM);
  else if (str_compare(flag,"tetragonal"))
    return(CELL_TYPE_TETR);
  else if (str_compare(flag,"hexagonal"))
    return(CELL_TYPE_HEXA);
  else if (str_compare(flag,"cubic"))
    return(CELL_TYPE_CUBC);
  else
    msg_error_f("unknown cell type \"%s\"",1,flag);
  return(0);
  }

/* convert internal cell type ID to corresponding flag
   
   t    - the internal ID
   flag - storage for the flag */
char* cell_type_name(short t, char *flag) {
  flag[0] = '\0';
  switch (t) {
    case CELL_TYPE_TRIC: sprintf(flag,"triclinic");    break;
    case CELL_TYPE_MONO: sprintf(flag,"monoclinic");   break;
    case CELL_TYPE_ORTH: sprintf(flag,"orthorhombic"); break;
    case CELL_TYPE_RHOM: sprintf(flag,"rhombohedral"); break;
    case CELL_TYPE_TETR: sprintf(flag,"tetragonal");   break;
    case CELL_TYPE_HEXA: sprintf(flag,"hexagonal");    break;
    case CELL_TYPE_CUBC: sprintf(flag,"cubic");        break;
    }
  return(flag);
  }

/* deduce and set cell type from cell vectors
 
   c - the cell data */
void cell_type_set(struct cell *c) {
  /* orthogonal systems */
  if (fabs(c->angle[0]-90.0)<CELL_CHECK_ACC &&
      fabs(c->angle[1]-90.0)<CELL_CHECK_ACC &&
      fabs(c->angle[2]-90.0)<CELL_CHECK_ACC) {
    /* cubic cell */
    if (fabs(c->side[0]-c->side[1])<CELL_CHECK_ACC &&
        fabs(c->side[0]-c->side[2])<CELL_CHECK_ACC &&
        fabs(c->side[1]-c->side[2])<CELL_CHECK_ACC)
      c->type = CELL_TYPE_CUBC;
    /* tetragonal cell */
    else if (fabs(c->side[0]-c->side[1])<CELL_CHECK_ACC ||
             fabs(c->side[0]-c->side[2])<CELL_CHECK_ACC ||
             fabs(c->side[1]-c->side[2])<CELL_CHECK_ACC)
      c->type = CELL_TYPE_TETR;
    /* orthorhombic */
    else
      c->type = CELL_TYPE_ORTH;
    }
  /* hexagonal cell */
  else if ((fabs(c->side[0]-c->side[1])<CELL_CHECK_ACC ||
            fabs(c->side[0]-c->side[2])<CELL_CHECK_ACC ||
            fabs(c->side[1]-c->side[2])<CELL_CHECK_ACC) &&
           (fabs(c->angle[0]-120.0)<CELL_CHECK_ACC ||
            fabs(c->angle[1]-120.0)<CELL_CHECK_ACC ||
            fabs(c->angle[2]-120.0)<CELL_CHECK_ACC))
    c->type = CELL_TYPE_HEXA;
  /* rhombohedral cell */
  else if (fabs(c->angle[0]-c->angle[1])<CELL_CHECK_ACC &&
           fabs(c->angle[0]-c->angle[2])<CELL_CHECK_ACC &&
           fabs(c->angle[1]-c->angle[2])<CELL_CHECK_ACC &&
           fabs(c->side[0]-c->side[1])<CELL_CHECK_ACC &&
           fabs(c->side[0]-c->side[2])<CELL_CHECK_ACC &&
           fabs(c->side[1]-c->side[2])<CELL_CHECK_ACC)
    c->type = CELL_TYPE_RHOM;
  /* monoclinic cell */
  else if ((fabs(c->angle[0]-90.0)<CELL_CHECK_ACC &&
            fabs(c->angle[1]-90.0)<CELL_CHECK_ACC) ||
           (fabs(c->angle[0]-90.0)<CELL_CHECK_ACC &&
            fabs(c->angle[2]-90.0)<CELL_CHECK_ACC) ||
           (fabs(c->angle[1]-90.0)<CELL_CHECK_ACC &&
            fabs(c->angle[2]-90.0)<CELL_CHECK_ACC))
    c->type = CELL_TYPE_MONO;
  /* triclinic cell */
  else
    c->type = CELL_TYPE_TRIC;
  /* set transformation matrices */
  cell_tmat_set(c);
  }

/* -------------------------------------------------------------------------- */
