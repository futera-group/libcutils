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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* write list of atoms to LAMMPS data file
 
   a - atomic data
   t - atom style
   n - number of atoms
   f - open file stream */
void lammps_atom_write(struct lammps_atom *a, short t, unsigned n, FILE *f) {
  unsigned i;
  if (a && n) {
    fprintf(f,"\nAtoms # %s\n\n",lammps_atom_type_sym(t));
    switch (t) {
      case LAMMPS_ATOM_CHRG:
        for (i=0; i<n; i++)
          fprintf(f,"%10d%5d%10.4f%14.6e%14.6e%14.6e\n",
            i+1,a[i].type+1,a[i].charge,a[i].crd[0],a[i].crd[1],a[i].crd[2]);
        break;
      case LAMMPS_ATOM_FULL:
        for (i=0; i<n; i++)
          fprintf(f,"%10d%10d%5d%10.4f%14.6e%14.6e%14.6e\n",
            i+1,a[i].mol+1,a[i].type+1,a[i].charge,
            a[i].crd[0],a[i].crd[1],a[i].crd[2]);
        break;
      default:
        msg_error_f("output of atom style (%d) is not supported",1,t);
      }
    }
  }

/* write list of atomic velocities to LAMMPS data file
 
   a - atomic data
   t - atom style
   n - number of atoms
   f - open file stream */
void lammps_atom_write_vel(struct lammps_atom *a, short t,
  unsigned n, FILE *f) {
  unsigned i;
  if (a && n && a[0].vel) {
    fprintf(f,"\nVelocities\n\n");
    for (i=0; i<n; i++) {
      fprintf(f,"%10d",i+1);
      fprintf(f,"%14.6e%14.6e%14.6e",
        a[i].vel[0],a[i].vel[1],a[i].vel[2]);
      switch (t) {
        case LAMMPS_ATOM_ELEC:
          fprintf(f,"%14.6e",a[i].rvl[0]);
          break;
        case LAMMPS_ATOM_ELLP:
          fprintf(f,"%14.6e%14.6e%14.6e",
            a[i].ang[0],a[i].ang[1],a[i].ang[2]);
          break;
        case LAMMPS_ATOM_SPHR:
          fprintf(f,"%14.6e%14.6e%14.6e",
            a[i].avl[0],a[i].avl[1],a[i].avl[2]);
          break;
        }
      fprintf(f,"\n");
      }
    }
  }

/* -------------------------------------------------------------------------- */
