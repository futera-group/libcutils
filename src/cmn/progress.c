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
#include <stdlib.h>
#include "cmn/message.h"
#include "cmn/progress.h"

/* -------------------------------------------------------------------------- */

/* allocate new progress bar */
struct progress *progress_new(void) {
  struct progress *p = NULL;
  /* allocate memory */
  p = (struct progress*)malloc(sizeof(struct progress));
  if (!p)
    msg_error("cannot allocate memory for new progress bar",1);
  /* initialization */
  p->step_id = 0;
  p->p_break = 2;
  p->p_step = 2;
  p->sign_reg = '-';
  p->sign_dec = '*';
  p->n_dt_total = 0;
  p->n_ps_printed = 0;
  p->n_ps_total = 50;
  p->flush = 0;
  p->file = stdout;
  return(p);
  }

/* clean memory allocated for progress bar
 
   p - pointer to progress bar data struct */
void progress_free(struct progress *p) {
  if (p)
    free(p);
  }

/* initialization progress bar print-out
 
   p - pointer to progress bar data struct */
void progress_init(struct progress *p, long unsigned n, FILE *f) {
  p->step_id = 0;
  p->n_dt_total = n;
  p->file = (f ? f : stdout);
  p->flush = ((p->file==stdout || p->file==stderr) ? 1 : 0);
  p->p_step = 100/p->n_ps_total;
  p->p_break = p->p_step;
  fprintf(p->file," |");
  if (p->flush)
    fflush(p->file);
  }

/* finish the progress bar print-out 
 
   p - pointer to progress data struct */
void progress_finish(struct progress *p) {
  unsigned i;
  for (i=p->n_ps_printed; i<(p->n_ps_total-1); i++)
    fprintf(p->file,"%c",(!((i+1)%5) ? '*' : '-'));
  fprintf(p->file,"|\n");
  fflush(p->file);
  }

/* calculate progress percentage and increase progress bar if needed
 
   p  - pointer to progress bar data struct
   id - data step ID */
void progress_step(struct progress *p, long unsigned id) {
  unsigned pg;
  pg = round((100.0*id)/((double)p->n_dt_total));
  while (pg<100 && pg>=p->p_break) {
    fprintf(p->file,"%c",((p->p_break)%10==0 ? '*' : '-'));
    if (p->flush)
      fflush(p->file);
    p->p_break += p->p_step;
    p->n_ps_printed++;
    }
  p->step_id = id;
  }

/* add one step and increase the progress bar if needed
 
   p  - pointer to progress bar data struct */
void progress_step_one(struct progress *p) {
  p->step_id++;
  progress_step(p,p->step_id);
  }

/* -------------------------------------------------------------------------- */
