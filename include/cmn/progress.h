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

#ifndef ZF_LIB_CMN_PROGRESS_H
#define ZF_LIB_CMN_PROGRESS_H

#include <stdio.h>

/* -------------------------------------------------------------------------- */

/* progress bar data struct */
struct progress {
  long unsigned step_id;     /* current step */
  unsigned p_break;          /* progress bar break point */
  unsigned p_step;           /* progress bar step */
  char sign_reg;             /* regular sign */
  char sign_dec;             /* decade sign */
  long unsigned n_dt_total;  /* total number of data values */
  unsigned n_ps_printed;     /* number of printed progress signs */
  unsigned n_ps_total;       /* total number of progress signs */
  short flush;               /* flush stream indicator */
  FILE *file;                /* output stream */
  };

/* -------------------------------------------------------------------------- */

/* allocate new progress bar */
struct progress *progress_new(void);
/* clean memory allocated for progress bar */
void progress_free(struct progress*);

/* initialization progress bar print-out */
void progress_init(struct progress*, long unsigned, FILE*);
/* finish the progress bar print-out */
void progress_finish(struct progress*);

/* calculate progress percentage and increase progress bar if needed */
void progress_step(struct progress*, long unsigned);
/* add one step and increase the progress bar if needed */
void progress_step_one(struct progress*);

/* -------------------------------------------------------------------------- */

#endif
