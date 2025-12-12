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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* write structure and force-field data to LAMMPS data file
 
   d    - the structure data storage
   name - name of the file */
void lammps_dat_write(struct lammps_dat *d, char *name) {
  FILE *f = stdout;
  /* open the file */
  if (name && name[0])
    f = file_open(name,"w");
  /* title */
  fprintf(f,"%s\n",d->header);
  fprintf(f,"\n");
  /* number of particles / types */
  if (d->n_atoms)
    fprintf(f,"%d atoms\n",d->n_atoms);
  if (d->n_atom_types)
    fprintf(f,"%d atom types\n",d->n_atom_types);
  if (d->n_bonds)
    fprintf(f,"%d bonds\n",d->n_bonds);
  if (d->n_bond_types)
    fprintf(f,"%d bond types\n",d->n_bond_types);
  if (d->n_angles)
    fprintf(f,"%d angles\n",d->n_angles);
  if (d->n_angle_types)
    fprintf(f,"%d angle types\n",d->n_angle_types);
  if (d->n_dihedrals)
    fprintf(f,"%d dihedrals\n",d->n_dihedrals);
  if (d->n_dihed_types)
    fprintf(f,"%d dihedral types\n",d->n_dihed_types);
  if (d->n_impropers)
    fprintf(f,"%d impropers\n",d->n_impropers);
  if (d->n_improp_types)
    fprintf(f,"%d improper types\n",d->n_improp_types);
  if (d->n_extra_bonds)
    fprintf(f,"%d extra bond per atom\n",d->n_extra_bonds);
  if (d->n_extra_angles)
    fprintf(f,"%d extra angle per atom\n",d->n_extra_angles);
  if (d->n_extra_dihedrals)
    fprintf(f,"%d extra dihedral per atom\n",d->n_extra_dihedrals);
  if (d->n_extra_impropers)
    fprintf(f,"%d extra improper per atom\n",d->n_extra_impropers);
  if (d->n_extra_specials)
    fprintf(f,"%d extra special per atom\n",d->n_extra_specials);
  if (d->n_ellipsoids)
    fprintf(f,"%d ellipsoids\n",d->n_ellipsoids);
  if (d->n_lines)
    fprintf(f,"%d lines\n",d->n_lines);
  if (d->n_triangles)
    fprintf(f,"%d triangles\n",d->n_triangles);
  if (d->n_bodies)
    fprintf(f,"%d bodies\n",d->n_bodies);
  fprintf(f,"\n");
  /* simulation box */
  lammps_box_write(d->box,f);
  /* particle types */
  lammps_type_write(d->type,d->n_atom_types,f);
  /* potentials */
  lammps_pot_write(d->pot_bond,LAMMPS_POT_BOND,1,d->n_bond_types,f);
  lammps_pot_write(d->pot_angle,LAMMPS_POT_ANGL,1,d->n_angle_types,f);
  lammps_pot_write(d->pot_dihed,LAMMPS_POT_DIHE,1,d->n_dihed_types,f);
  /* atoms & velocities */
  lammps_atom_write(d->atom,d->atom_style,d->n_atoms,f);
  lammps_atom_write_vel(d->atom,d->atom_style,d->n_atoms,f);
  /* bonding-potential IDs */
  lammps_pot_write(d->bond,LAMMPS_POT_BOND,0,d->n_bonds,f);
  lammps_pot_write(d->angle,LAMMPS_POT_ANGL,0,d->n_angles,f);
  lammps_pot_write(d->dihed,LAMMPS_POT_DIHE,0,d->n_dihedrals,f);
  lammps_pot_write(d->impr,LAMMPS_POT_IMPR,0,d->n_impropers,f);
  /* pair potentials */
  lammps_pot_write(d->pot_pair,LAMMPS_POT_PAIR,0,d->n_atom_types,f);
  /* close the file */
  if (name && name[0]) 
    file_close(f);
  }

/* -------------------------------------------------------------------------- */
