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

#include <string.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read atomic data from cp2k output file
 
   d - pointer to cp2k data struct
   f - open file stream */
void cp2k_log_read_atoms(struct cp2k_dat *d, FILE *f) {
  unsigned i;
  char *line;
  /* find the atom data block */
  rewind(f);
  line = str_ffind_new(f," ATOMIC COORDINATES IN");
  if (line) {
    line = str_free(line);
    /* skip header */
    for (i=0; i<3; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      if (strstr(line,"Atom")) {
        line = str_free(line);
        break;
        }
      line = str_free(line);
      }
    /* read atomic data */
    if (!d->atom)
      d->atom = cp2k_atom_new(d->n_atoms);
    for (i=0; i<d->n_atoms; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      /* skip blank line present in output of older CP2K versions */
      if (i==0 && str_length(line) <= 1) {
        line = str_free(line);
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of cp2k log file",1);
        }
      /* atomic data */
      if (sscanf(line,"%*d%d%*s%d%lf%lf%lf%lf%lf",
          &(d->atom[i].kind),&(d->atom[i].num),
          &(d->atom[i].crd[0]),&(d->atom[i].crd[1]),&(d->atom[i].crd[2]),
          &(d->atom[i].charge[CHRG_EFF]),&(d->atom[i].mass))!=7)
        msg_error("invalid atom data",1);
      if (d->atom[i].kind<1 || d->atom[i].kind>d->n_atom_kinds)
        msg_error_f("invalid kind ID (%d) of atom #%d",1,d->atom[i].kind,i+1);
      d->atom[i].kind--;
      line = str_free(line);
      }
    }
  }

/* -------------------------------------------------------------------------- */
