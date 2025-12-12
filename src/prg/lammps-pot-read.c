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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* read potential specification from LAMMPS data file 
 
   pot_type  - type of the potential
   pot_style - potential style
   read_id   - read type ID from file
   n_ids     - number of atom IDs
   n_coeffs  - number of coefficients
   n_pots    - number of potential data lines 
   f         - open file stream */
struct lammps_pot* lammps_pot_read(short pot_type, short pot_style,
  short read_id, short read_tp, unsigned n_ids, unsigned n_coeffs,
  unsigned n_pots, FILE *f) {
  char *line,**s;
  unsigned i,ns,nt,pot_id,ir,id = 0;
  short *check;
  struct lammps_pot *p;
  /* initialization */
  nt = n_ids+n_coeffs+(read_tp ? 1 : 0);
  p = lammps_pot_new(n_pots);
  if (read_id) {
    check = vec_sialloc(n_pots);
    vec_siset(check,0,n_pots);
    nt++;
    }
  /* read data from file */
  for (line=str_read_line_new(f); id<n_pots && line;
       line=str_free(line),line=str_read_line_new(f)) {
    str_trim(line);
    /* blank line */
    if (str_length(line)<1)
      continue;
    /* check data */
    s = str_split(line,' ',&ns);
    if (ns!=nt)
      msg_error_f("inconsistent number of data fields (%d/%d) "
        "in potential type %d(%d)",1,ns,nt,pot_type,pot_style);
    ir = 0;
    /* potential ID */
    if (read_id) {
      if (sscanf(s[ir++],"%u",&pot_id)!=1)
        msg_error_f("invalid format of potential ID in type %d(%d)",1,
          pot_type,pot_style);
      if (pot_id<1 || pot_id>n_pots)
        msg_error_f("invalid value of potential ID in type %d(%d)",1,
          pot_type,pot_style);
      pot_id--;
      }
    else
      pot_id = id;
    /* save data */
    p[pot_id].style = pot_style;
    if (read_tp) {
      if (sscanf(s[ir++],"%u",&(p[pot_id].type))!=1)
        msg_error_f("invalid format of type ID in potential type %d(%d)",1,
          pot_type,pot_style);
      if (p[pot_id].type<1)
        msg_error_f("invalid value of type ID in potential type %d(%d)",1,
          pot_type,pot_style);
      p[pot_id].type--;
      }
    p[pot_id].n_ids = n_ids;
    p[pot_id].id = vec_ualloc(p[pot_id].n_ids);
    for (i=0; i<p[pot_id].n_ids; i++)
      if (sscanf(s[ir++],"%u",&(p[pot_id].id[i]))!=1)
        msg_error_f("invalid format of atom ID in potential type %d(%d)",1,
          pot_type,pot_style);
    p[pot_id].n_coeffs = n_coeffs;
    p[pot_id].coeff = vec_falloc(p[pot_id].n_coeffs);
    for (i=0; i<p[pot_id].n_coeffs; i++)
      if (sscanf(s[ir++],"%lf",&(p[pot_id].coeff[i]))!=1)
        msg_error_f("invalid format of coefficient in potential type %d(%d)",1,
          pot_type,pot_style);
    /* next potential */
    s = vec_sfree(s,ns);
    if (read_id) 
      check[pot_id] = 1;
    id++;
    }
  /* sanity check */
  if (id<n_pots)
    msg_error_f("unexpected error while reading specification"
      " of potential type %d(%d)",1,pot_type,pot_style);
  if (read_id) {
    for (i=0; i<n_pots; i++)
      if (!check[i])
        msg_error_f("potential #%d of %d(%d) not assigned",1,
          i+1,pot_type,pot_style);
    vec_sifree(check);
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
