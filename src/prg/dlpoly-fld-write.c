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
#include <cmn/file.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* write force-field parameters to external file
 
   d    - the force-field data
   name - name of the file */
void dlpoly_fld_write(struct dlpoly_fld *d, char *name) {
  unsigned i,j,n;
  FILE *f = stdout;
  if (name) 
    f = file_open(name,"w");
  /* header */
  fprintf(f,"%s\n",d->header);
  /* units */
  fprintf(f,"UNITS %s\n",dlpoly_fld_unit_sym(d->units));
  if (name)
    file_close(f);
  /* molecules */
  fprintf(f,"MOLECULES %d\n",d->n_molecules);
  for (i=0; i<d->n_molecules; i++) {
    fprintf(f,"%s\n",d->mol[i].name);
    fprintf(f,"NUMMOLS %d\n",d->mol[i].n_molecules);
    /* atoms */
    fprintf(f,"ATOMS %d\n",d->mol[i].n_atoms);
    for (j=0, n=0; j<d->mol[i].n_atoms && n<d->mol[i].n_atoms; j++) {
      fprintf(f,"  %-8s %8.3f %10.6f %2d %2d\n",
        d->mol[i].atom[j].name,
        d->mol[i].atom[j].mass,
        d->mol[i].atom[j].charge,
        d->mol[i].atom[j].n_atoms,
        d->mol[i].atom[j].frozen);
      n = n + d->mol[i].atom[j].n_atoms;
      }
    /* bonds */
    if (d->mol[i].n_bonds) {
      fprintf(f,"BONDS %d\n",d->mol[i].n_bonds);
      for (j=0; j<d->mol[i].n_bonds; j++) {
        n = dlpoly_prm_bond_pot_nvar(d->mol[i].bond[j].pot);
        fprintf(f,"  %-8s %3d %3d",
          dlpoly_prm_bond_pot_sym(d->mol[i].bond[j].pot,0),
          d->mol[i].bond[j].id1+1,
          d->mol[i].bond[j].id2+1);
        if (n>0) fprintf(f," %8.3f",d->mol[i].bond[j].v1);
        if (n>1) fprintf(f," %8.3f",d->mol[i].bond[j].v2);
        if (n>2) fprintf(f," %8.3f",d->mol[i].bond[j].v3);
        if (n>3) fprintf(f," %8.3f",d->mol[i].bond[j].v4);
        fprintf(f,"\n");
        }
      }
    /* angles */
    if (d->mol[i].n_angles) {
      fprintf(f,"ANGLES %d\n",d->mol[i].n_angles);
      for (j=0; j<d->mol[i].n_angles; j++) {
        n = dlpoly_prm_angle_pot_nvar(d->mol[i].angle[j].pot);
        fprintf(f,"  %-8s %3d %3d %3d",
          dlpoly_prm_angle_pot_sym(d->mol[i].angle[j].pot,0),
          d->mol[i].angle[j].id1+1,
          d->mol[i].angle[j].id2+1,
          d->mol[i].angle[j].id3+1);
        if (n>0) fprintf(f," %8.3f",d->mol[i].angle[j].v1);
        if (n>1) fprintf(f," %8.3f",d->mol[i].angle[j].v2);
        if (n>2) fprintf(f," %8.3f",d->mol[i].angle[j].v3);
        if (n>3) fprintf(f," %8.3f",d->mol[i].angle[j].v4);
        fprintf(f,"\n");
        }
      }
    /* dihedral angles */
    if (d->mol[i].n_dihedrals) {
      fprintf(f,"DIHEDRALS %d\n",d->mol[i].n_dihedrals);
      for (j=0; j<d->mol[i].n_dihedrals; j++) {
        n = dlpoly_prm_dihed_pot_nvar(d->mol[i].dihed[j].pot);
        fprintf(f,"  %-8s %3d %3d %3d %3d",
          dlpoly_prm_dihed_pot_sym(d->mol[i].dihed[j].pot),
          d->mol[i].dihed[j].id1+1,
          d->mol[i].dihed[j].id2+1,
          d->mol[i].dihed[j].id3+1,
          d->mol[i].dihed[j].id4+1);
        if (n>0) fprintf(f," %8.3f",d->mol[i].dihed[j].v1);
        if (n>1) fprintf(f," %8.3f",d->mol[i].dihed[j].v2);
        if (n>2) fprintf(f," %8.3f",d->mol[i].dihed[j].v3);
        if (n>3) fprintf(f," %8.3f",d->mol[i].dihed[j].v4);
        if (n>3) fprintf(f," %8.3f",d->mol[i].dihed[j].v5);
        if (n>3) fprintf(f," %8.3f",d->mol[i].dihed[j].v6);
        if (n>4) fprintf(f," %8.3f",d->mol[i].dihed[j].v7);
        fprintf(f,"\n");
        }
      }
    fprintf(f,"FINISH\n");
    }
  /* VDW potentials */
  if (d->n_vdws) {
    fprintf(f,"VDW %d\n",d->n_vdws);
    for (i=0; i<d->n_vdws; i++) {
      n = dlpoly_prm_vdw_pot_nvar(d->vdw[i].pot);
      fprintf(f,"  %-8s %-8s %-8s",
        d->vdw[i].at1,d->vdw[i].at2,dlpoly_prm_vdw_pot_sym(d->vdw[i].pot));
      if (n>0) fprintf(f," %12.3f",d->vdw[i].v1);
      if (n>1) fprintf(f," %12.3f",d->vdw[i].v2);
      if (n>2) fprintf(f," %12.3f",d->vdw[i].v3);
      if (n>3) fprintf(f," %12.3f",d->vdw[i].v4);
      if (n>4) fprintf(f," %12.3f",d->vdw[i].v5);
      fprintf(f,"\n");
      }
    }
  fprintf(f,"CLOSE\n");
  }

/* -------------------------------------------------------------------------- */
