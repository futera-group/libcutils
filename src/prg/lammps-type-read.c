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

#include <stdio.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* read atomic types from LAMMPS data file
 
   n - number of atomic types
   f - open file stream */
struct lammps_type *lammps_type_read(unsigned n, FILE *f) {
  char *line,**s;
  unsigned i,ns,tp_id,id = 0;
  short *check;
  struct lammps_type *m;
  /* initalization */
  m = lammps_type_new(n);
  check = vec_sialloc(n);
  vec_siset(check,0,n);
  /* read data from file */
  for (line=str_read_line_new(f); id<n && line;
       line=str_free(line),line=str_read_line_new(f)) {
    str_trim(line);
    /* blank line */
    if (str_length(line)<1)
      continue;
    /* check data */
    s = str_split(line,' ',&ns);
    if (ns<2)
      msg_error_f("inconsistent number of atom type coefficients (%d/2)",1,ns);
    /* type ID */
    if (sscanf(s[0],"%u",&tp_id)!=1)
      msg_error("invalid format of atomic type ID",1);
    if (tp_id<1 || tp_id>n)
      msg_error_f("invalid value of atomic type ID (%d)",1,tp_id);
    tp_id--;
    /* atomic mass */
    if (sscanf(s[1],"%lf",&(m[tp_id].mass))!=1)
      msg_error_f("invalid format of atomic mass in type #%d",1,tp_id+1);
    /* atomic number & name */
    m[tp_id].num = atom_num_mass(m[tp_id].mass);
    if (ns>3 && str_compare(s[2],"#"))
      m[tp_id].name = str_copy_new(s[3]);
    else
      m[tp_id].name = str_copy_new(atom_name(m[tp_id].num));
    /* next type */
    s = vec_sfree(s,ns);
    check[tp_id] = 1;
    id++;
    }
  /* sanity check */
  if (id<n)
    msg_error("unexpected error while reading atom types",1);
  for (i=0; i<n; i++)
    if (!check[i])
      msg_error_f("atom type #%d not assigned",1,i+1);
  vec_sifree(check);
  return(m);
  }

/* -------------------------------------------------------------------------- */
