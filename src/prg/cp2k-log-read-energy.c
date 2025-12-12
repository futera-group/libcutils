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

#include <cmn/message.h>
#include <cmn/stack.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read energies from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_energy(struct cp2k_dat *d, FILE *f) {
  char *line;
  double en;
  struct stack *t;
  /* initialization */
  t = stack_alloc();
  rewind(f);
  /* read all Fermi-energy records */
  for (line=str_ffind_b_new(f,"  Fermi energy:"); line;
       line=str_free(line),line=str_ffind_b_new(f,"  Fermi energy:")) {
    if (str_length(line)<31 || sscanf(line+30,"%lf",&en)!=1)
      msg_error("cannot read fermi energy\n",1);
    stack_fadd(t,en);
    }
  /* save final Fermi levels */
  if (t->num)
    stack_fget(t,&(d->fermi));
  /* clean memory */
  stack_free(t);
  }

/* -------------------------------------------------------------------------- */
