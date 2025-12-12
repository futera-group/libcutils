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
#include <cmn/message.h>
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for gaussian log-file basis-set data structure */
struct gauss_log_bs* gauss_log_bs_new(void) {
  unsigned i;
  struct gauss_log_bs *d = NULL;
  /* memory allocation */
  d = (struct gauss_log_bs*)malloc(sizeof(struct gauss_log_bs));
  if (!d)
    msg_error("cannot allocate memory for gaussian log-file basis-set data",1);
  /* initialization */
  d->at_id = 0;
  d->at_num = 0;
  d->shell = basis_shell_new(1);
  for (i=0; i<3; i++)
    d->coord[i] = 0.0;
  return(d);
  }

/* free memory allocate for log-file basis-set data structure

   d - the basis-set data structure */
struct gauss_log_bs* gauss_log_bs_free(struct gauss_log_bs *d) {
  if (d) {
    d->shell = basis_shell_free(d->shell,1);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
