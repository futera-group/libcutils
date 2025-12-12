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

/* read atomic velocities from LAMMPS data file

   t  - atom style
   nv - number of expected coefficients 
   n  - number of atoms
   f  - open file stream */
void lammps_atom_read_vel(struct lammps_atom *a, short t, unsigned nv, 
  unsigned n, FILE *f) {
  char *line,**s;
  unsigned i,ns,at_id,id = 0;
  short *check;
  /* initalization */
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
    if (ns<nv)
      msg_error_f("inconsistent number of velocity coefficients (%d/%d) "
        "in atom style %d",1,ns,nv,t);
    /* atom ID */
    if (sscanf(s[0],"%u",&at_id)!=1)
      msg_error("invalid format of velocity atom ID",1);
    if (at_id<1 || at_id>n)
      msg_error_f("invalid value of velocity atom ID (%d)",1,at_id);
    at_id--;
    /* velocities */
    a[at_id].vel = vec_falloc(3);
    if (sscanf(s[1],"%lf",&(a[at_id].vel[0]))!=1)
      msg_error_f("invalid format of atom #%d X velocity",1,at_id+1);
    if (sscanf(s[2],"%lf",&(a[at_id].vel[1]))!=1)
      msg_error_f("invalid format of atom #%d Y velocity",1,at_id+1);
    if (sscanf(s[3],"%lf",&(a[at_id].vel[2]))!=1)
      msg_error_f("invalid format of atom #%d Z velocity",1,at_id+1);
    /* additional data */
    switch (t) {
      case LAMMPS_ATOM_ELEC:
        a[at_id].rvl = vec_falloc(1);
        if (sscanf(s[4],"%lf",&(a[at_id].rvl[0]))!=1)
          msg_error_f("invalid format of atom #%d radial velocity",1,at_id+1);
        break;
      case LAMMPS_ATOM_ELLP:
        a[at_id].ang = vec_falloc(3);
        if (sscanf(s[4],"%lf",&(a[at_id].ang[0]))!=1)
          msg_error_f("invalid format of atom #%d X ang. momentum",1,at_id+1);
        if (sscanf(s[5],"%lf",&(a[at_id].ang[1]))!=1)
          msg_error_f("invalid format of atom #%d Y ang. momentum",1,at_id+1);
        if (sscanf(s[6],"%lf",&(a[at_id].ang[2]))!=1)
          msg_error_f("invalid format of atom #%d Z ang. momentum",1,at_id+1);
        break;
      case LAMMPS_ATOM_SPHR:
        a[at_id].avl = vec_falloc(3);
        if (sscanf(s[4],"%lf",&(a[at_id].avl[0]))!=1)
          msg_error_f("invalid format of atom #%d X ang. velocity",1,at_id+1);
        if (sscanf(s[5],"%lf",&(a[at_id].avl[1]))!=1)
          msg_error_f("invalid format of atom #%d X ang. velocity",1,at_id+1);
        if (sscanf(s[6],"%lf",&(a[at_id].avl[2]))!=1)
          msg_error_f("invalid format of atom #%d X ang. velocity",1,at_id+1);
        break;
      }
    /* next atom */
    s = vec_sfree(s,ns);
    check[at_id] = 1;
    id++;
    }
  /* sanity check */
  if (id<n)
    msg_error("unexpected error while reading atom velocities",1);
  for (i=0; i<n; i++)
    if (!check[i])
      msg_error_f("atom velocity #%d not assigned",1,i+1);
  vec_sifree(check);
  }

/* read list of atoms from LAMMPS data file
 
   t  - atom style
   nv - number of expected coefficients 
   n  - number of atoms
   f  - open file stream */
struct lammps_atom* lammps_atom_read(short t, unsigned nv,
  unsigned n, FILE *f) {
  char *line,**s;
  unsigned i,ns,at_id,id = 0;
  short *check;
  struct lammps_atom *a;
  /* initalization */
  a = lammps_atom_new(n);
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
    if (ns<nv)
      msg_error_f("inconsistent number of coefficients (%d/%d) "
        "in atom style %d",1,ns,nv,t);
    /* atom ID */
    if (sscanf(s[0],"%u",&at_id)!=1)
      msg_error("invalid format of atom ID",1);
    if (at_id<1 || at_id>n)
      msg_error_f("invalid value of atom ID (%d)",1,at_id);
    at_id--;
    /* memory allocation */
    a[at_id].crd = vec_falloc(3);
    /* style-specific data assignment */
    switch (t) {
      case LAMMPS_ATOM_CHRG:
        if (sscanf(s[1],"%u",&(a[at_id].type))!=1)
          msg_error_f("invalid format of atom #%d type",1,at_id+1);
        if (a[at_id].type<1)
          msg_error_f("invalid value of  atom #%d type",1,at_id+1);
        a[at_id].type--;
        if (sscanf(s[2],"%lf",&(a[at_id].charge))!=1)
          msg_error_f("invalid format of atom #%d charge",1,at_id+1);
        if (sscanf(s[3],"%lf",&(a[at_id].crd[0]))!=1)
          msg_error_f("invalid format of atom #%d X coordinate",1,at_id+1);
        if (sscanf(s[4],"%lf",&(a[at_id].crd[1]))!=1)
          msg_error_f("invalid format of atom #%d Y coordinate",1,at_id+1);
        if (sscanf(s[5],"%lf",&(a[at_id].crd[2]))!=1)
          msg_error_f("invalid format of atom #%d Z coordinate",1,at_id+1);
        break;
      case LAMMPS_ATOM_FULL:
        if (sscanf(s[1],"%d",&(a[at_id].mol))!=1)
          msg_error_f("invalid format of molecular ID of atom #%d",1,at_id+1);
        if (a[at_id].mol<1)
          msg_error_f("invalid value of molecular ID of atom #%d",1,at_id+1);
        a[at_id].mol--;
        if (sscanf(s[2],"%u",&(a[at_id].type))!=1)
          msg_error_f("invalid format of atom #%d type",1,at_id+1);
        if (a[at_id].type<1)
          msg_error_f("invalid value of  atom #%d type",1,at_id+1);
        a[at_id].type--;
        if (sscanf(s[3],"%lf",&(a[at_id].charge))!=1)
          msg_error_f("invalid format of atom #%d charge",1,at_id+1);
        if (sscanf(s[4],"%lf",&(a[at_id].crd[0]))!=1)
          msg_error_f("invalid format of atom #%d X coordinate",1,at_id+1);
        if (sscanf(s[5],"%lf",&(a[at_id].crd[1]))!=1)
          msg_error_f("invalid format of atom #%d Y coordinate",1,at_id+1);
        if (sscanf(s[6],"%lf",&(a[at_id].crd[2]))!=1)
          msg_error_f("invalid format of atom #%d Z coordinate",1,at_id+1);
        break;
      default:
        msg_error_f("atom style %d not supported",1,t);
      }
    /* next atom */
    s = vec_sfree(s,ns);
    check[at_id] = 1;
    id++;
    }
  /* sanity check */
  if (id<n)
    msg_error("unexpected error while reading atoms",1);
  for (i=0; i<n; i++)
    if (!check[i])
      msg_error_f("atom #%d not assigned",1,i+1);
  vec_sifree(check);
  return(a);
  }

/* -------------------------------------------------------------------------- */
