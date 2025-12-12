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

/* read Mulliken charges from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_charges_mulliken(struct cp2k_dat *d, FILE *f) {
  unsigned i;
  short spin;
  char *line = NULL;
  /* Mulliken charge block */
  /* memory allocation */
  if (!d->atom)
    d->atom = cp2k_atom_new(d->n_atoms);
  /* header */
  for (i=0; i<2; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete mulliken-charge data block",1);
    }
  /* spin polarizaton */
  spin = 0;
  if (strstr(line,"Spin moment"))
    spin = 1;
  /* charges */
  for (i=0; i<d->n_atoms; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete mulliken-charge data block",1);
    if (spin) {
      if (sscanf(line,"%*d%*s%*d%*f%*f%lf",
          &(d->atom[i].charge[CHRG_MULL]))!=1)
        msg_error_f("cannot read mulliken charge of atom #%d",1,i+1);
      }
    else {
      if (sscanf(line,"%*d%*s%*d%*f%lf",
          &(d->atom[i].charge[CHRG_MULL]))!=1)
        msg_error_f("cannot read mulliken charge of atom #%d",1,i+1);
      }
    }
  line = str_free(line);
  }

/* read Lowdin charges from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_charges_lowdin(struct cp2k_dat *d, FILE *f) {
  int i;
  short spin;
  char *line = NULL;
  /* Lowdin charge block */
  /* memory allocation */
  if (!d->atom)
    d->atom = cp2k_atom_new(d->n_atoms);
  /* header */
  for (i=0; i<2; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete lowdin-charge data block",1);
    if (strstr(line,"WARNING in population_analyses"))
      i -= 4;
    }
  /* spin polarizaton */
  spin = 0;
  if (strstr(line,"Spin moment"))
    spin = 1;
  /* charges */
  for (i=0; i<d->n_atoms; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete lowdin-charge data block",1);
    if (spin) {
      if (sscanf(line,"%*d%*s%*d%*f%*f%lf",
          &(d->atom[i].charge[CHRG_LOWD]))!=1)
        msg_error_f("cannot read lowdin charge of atom #%d",1,i+1);
      }
    else {
      if (sscanf(line,"%*d%*s%*d%*f%lf",
          &(d->atom[i].charge[CHRG_LOWD]))!=1)
        msg_error_f("cannot read lowdin charge of atom #%d",1,i+1);
      }
    }
  line = str_free(line);
  }

/* read Hirshfeld charges from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_charges_hirshfeld(struct cp2k_dat *d, FILE *f) {
  unsigned i;
  short spin;
  char *line = NULL;
  /* Hirshfeld charge block */
  /* memory allocation */
  if (!d->atom)
    d->atom = cp2k_atom_new(d->n_atoms);
  /* header */
  for (i=0; i<2; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete hirshfeld-charge data block",1);
    }
  /* spin polarizaton */
  spin = 0;
  if (strstr(line,"Spin moment"))
    spin = 1;
  /* charges */
  for (i=0; i<d->n_atoms; i++) {
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("incomplete hirshfeld-charge data block",1);
    if (spin) {
      if (sscanf(line,"%*d%*s%*d%*f%*f%*f%*f%lf",
          &(d->atom[i].charge[CHRG_HIRSH]))!=1)
        msg_error_f("cannot read hirshfeld charge of atom #%d",1,i+1);
      }
    else {
      if (sscanf(line,"%*d%*s%*d%*f%*f%lf",
          &(d->atom[i].charge[CHRG_HIRSH]))!=1)
        msg_error_f("cannot read hirshfeld charge of atom #%d",1,i+1);
      }
    }
  line = str_free(line);
  }

/* read calculated atomic charges from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_charges(struct cp2k_dat *d, FILE *f) {
  char *line;
  /* initialization */
  rewind(f);
  /* search for charges */
  for (line=str_read_line_new(f); line;
       line=str_free(line), line=str_read_line_new(f)) {
    /* Mulliken charges */
    if (strstr(line,"Mulliken Population Analysis"))
      cp2k_log_read_charges_mulliken(d,f);
    /* Lowdin charges */
    if (strstr(line,"LOWDIN POPULATION ANALYSIS"))
      cp2k_log_read_charges_lowdin(d,f);
    /* Hirshfeld charges */
    if (strstr(line,"Hirshfeld Charges")) 
      cp2k_log_read_charges_hirshfeld(d,f);
    }
  }

/* -------------------------------------------------------------------------- */
