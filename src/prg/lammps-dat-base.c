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

#include <stdlib.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for LAMMPS structure and force-field data */
struct lammps_dat *lammps_dat_new(void) {
  struct lammps_dat *d = NULL;
  /* memory allocation */
  d = malloc(sizeof(struct lammps_dat));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS structure data",1);
  /* initialization */
  d->header = NULL;
  d->n_atoms = 0;
  d->n_bonds = 0;
  d->n_angles = 0;
  d->n_dihedrals = 0;
  d->n_impropers = 0;
  d->n_atom_types = 0;
  d->n_bond_types = 0;
  d->n_angle_types = 0;
  d->n_dihed_types = 0;
  d->n_improp_types = 0; 
  d->n_extra_bonds = 0;
  d->n_extra_angles = 0;
  d->n_extra_dihedrals = 0;
  d->n_extra_impropers = 0;
  d->n_extra_specials = 0;
  d->n_ellipsoids = 0;
  d->n_lines = 0;
  d->n_triangles = 0;
  d->n_bodies = 0;
  d->atom_style = 0;
  d->box = lammps_box_new();
  d->atom = NULL;
  d->bond = NULL;
  d->angle = NULL;
  d->dihed = NULL;
  d->impr = NULL;
  d->type = NULL;
  d->pot_bond = NULL;
  d->pot_angle = NULL;
  d->pot_dihed = NULL;
  d->pot_pair = NULL;
  d->trj = NULL;
  return(d);
  }

/* free memory allocated for LAMMPS structure and force-field data 
 
   d - the structure data */
struct lammps_dat *lammps_dat_free(struct lammps_dat *d) {
  if (d) {
    d->header = str_free(d->header);
    d->box = lammps_box_free(d->box);
    d->atom = lammps_atom_free(d->atom,d->n_atoms);
    d->bond = lammps_pot_free(d->bond,d->n_bonds);
    d->angle = lammps_pot_free(d->angle,d->n_angles);
    d->dihed = lammps_pot_free(d->dihed,d->n_dihedrals);
    d->impr = lammps_pot_free(d->impr,d->n_impropers);
    d->type = lammps_type_free(d->type,d->n_atom_types);
    d->pot_bond = lammps_pot_free(d->pot_bond,d->n_bond_types);
    d->pot_angle= lammps_pot_free(d->pot_angle,d->n_angle_types);
    d->pot_dihed= lammps_pot_free(d->pot_dihed,d->n_dihed_types);
    d->pot_pair= lammps_pot_free(d->pot_pair,d->n_atom_types);
    d->trj = lammps_trj_free(d->trj);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
