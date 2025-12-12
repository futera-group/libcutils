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

#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* initialize pointer struct according to pointer array

   p    - pointer to the amber pointer struct
   pvec - the pointer array */
void amber_pointer_init(struct amber_pointer *p, int *pvec) {
  p->natom = pvec[AMBER_POINTER_NATOM];
  p->ntype = pvec[AMBER_POINTER_NTYPES];
  p->natyp = pvec[AMBER_POINTER_NATYP];
  p->nres = pvec[AMBER_POINTER_NRES];
  p->nbond = pvec[AMBER_POINTER_NBOND];
  p->nangle = pvec[AMBER_POINTER_NANGLE];
  p->ndihed = pvec[AMBER_POINTER_NDIHED];
  p->nbondh = pvec[AMBER_POINTER_NBONH];
  p->nbonda = pvec[AMBER_POINTER_NBONA];
  p->nangleh = pvec[AMBER_POINTER_NTHETH];
  p->nanglea = pvec[AMBER_POINTER_NTHETA];
  p->ndihedh = pvec[AMBER_POINTER_NPHIH];
  p->ndiheda = pvec[AMBER_POINTER_NPHIA];
  p->next = pvec[AMBER_POINTER_NEXT];
  p->nhbond = pvec[AMBER_POINTER_NHBOND];
  }

/* -------------------------------------------------------------------------- */
