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
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* read structure and force-field data from LAMMPS data file
 
   d    - the structure data storage
   name - name of the file */
void lammps_dat_read(struct lammps_dat *d, char *name) {
  char *line,**v;
  unsigned n,nv;
  short tp,tc;
  FILE *f;
  /* open the file */
  f = file_open(name,"r");
  /* header */
  d->header = str_read_line_new(f);
  if (!d->header)
    msg_error_f("cannot read header line from \"%s\"",1,name);
  str_trim(d->header);
  /* keywords */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    str_trim(line);
    str_lowcase(line);
    v = str_split(line,' ',&n);
    /* blank line */
    if (n<1)
      continue;
    /* data lines */
    if (n>=6) {
      if (str_compare(v[3],"xy") && 
          str_compare(v[4],"xz") && 
          str_compare(v[5],"yz"))  {
        if (sscanf(v[0],"%lf",&(d->box->tilt[0]))!=1)
          msg_error("invalid value of XY tilting factor",1);
        if (sscanf(v[1],"%lf",&(d->box->tilt[1]))!=1)
          msg_error("invalid value of XZ tilting factor",1);
        if (sscanf(v[2],"%lf",&(d->box->tilt[2]))!=1)
          msg_error("invalid value of YZ tilting factor",1);
        continue;
        }
      }
    if (n>=5) {
      /* extra bonds */
      if (str_compare(v[1],"extra") && 
          str_compare(v[2],"bond") && 
          str_compare(v[3],"per") && 
          str_compare(v[4],"atom"))  {
        if (sscanf(v[0],"%u",&(d->n_extra_bonds))!=1)
          msg_error("invalid specification of number of extra bonds",1);
        continue;
        }
      /* extra angles */
      else if (str_compare(v[1],"extra") && 
          str_compare(v[2],"angle") && 
          str_compare(v[3],"per") && 
          str_compare(v[4],"atom")) {
        if (sscanf(v[0],"%u",&(d->n_extra_angles))!=1)
          msg_error("invalid specification of number of extra angles",1);
        continue;
        }
      /* extra dihedrals */
      else if (str_compare(v[1],"extra") && 
          str_compare(v[2],"dihedral") && 
          str_compare(v[3],"per") && 
          str_compare(v[4],"atom")) {
        if (sscanf(v[0],"%u",&(d->n_extra_dihedrals))!=1)
          msg_error("invalid specification of number of extra dihedrals",1);
        continue;
        }
      /* extra impropers */
      else if (str_compare(v[1],"extra") && 
          str_compare(v[2],"improper") && 
          str_compare(v[3],"per") && 
          str_compare(v[4],"atom")) {
        if (sscanf(v[0],"%u",&(d->n_extra_impropers))!=1)
          msg_error("invalid specification of number of extra impropers",1);
        continue;
        }
      /* extra special sites */
      else if (str_compare(v[1],"extra") && 
          str_compare(v[2],"special") && 
          str_compare(v[3],"per") && 
          str_compare(v[4],"atom")) {
        if (sscanf(v[0],"%u",&(d->n_extra_specials))!=1)
          msg_error("invalid specification of number of extra special sites",1);
        continue;
        }
      }
    if (n>=4) {
      /* box size in x-direction */
      if (str_compare(v[2],"xlo") && str_compare(v[3],"xhi")) {
        if (sscanf(v[0],"%lf",&(d->box->min[0]))!=1)
          msg_error("invalid value of PBC lower boundary in X direction",1);
        if (sscanf(v[1],"%lf",&(d->box->max[0]))!=1)
          msg_error("invalid value of PBC upper boundary in X direction",1);
        d->box->pbc[0] = 1;
        continue;
        }
      /* box size in y-direction */
      else if (str_compare(v[2],"ylo") && str_compare(v[3],"yhi")) {
        if (sscanf(v[0],"%lf",&(d->box->min[1]))!=1)
          msg_error("invalid value of PBC lower boundary in Y direction",1);
        if (sscanf(v[1],"%lf",&(d->box->max[1]))!=1)
          msg_error("invalid value of PBC upper boundary in Y direction",1);
        d->box->pbc[1] = 1;
        continue;
        }
      /* box size in z-direction */
      else if (str_compare(v[2],"zlo") && str_compare(v[3],"zhi")) {
        if (sscanf(v[0],"%lf",&(d->box->min[2]))!=1)
          msg_error("invalid value of PBC lower boundary in Z direction",1);
        if (sscanf(v[1],"%lf",&(d->box->max[2]))!=1)
          msg_error("invalid value of PBC upper boundary in Z direction",1);
        d->box->pbc[2] = 1;
        continue;
        }
      /* potential coefficients */
      else if (str_compare(v[1],"coeffs")) {
        tp = lammps_pot_type_id(v[0]);
        if (!tp)
          msg_error_f("unknown potential \"%s\"",1,v[0]);
        if (!str_compare(v[2],"#"))
          msg_error_f("style of \"%s\" potential is not specified",1,v[0]);
        switch (tp) {
          case LAMMPS_POT_BOND: 
            tc = lammps_pot_bond_id(v[3]);
            nv = lammps_pot_bond_nvar(tc);
            d->pot_bond = lammps_pot_read(tp,tc,1,0,0,nv,d->n_bond_types,f);
            break;
          case LAMMPS_POT_ANGL:
            tc = lammps_pot_angle_id(v[3]);
            nv = lammps_pot_angle_nvar(tc);
            d->pot_angle = lammps_pot_read(tp,tc,1,0,0,nv,d->n_angle_types,f);
            break;
          case LAMMPS_POT_DIHE:
            tc = lammps_pot_dihed_id(v[3]);
            nv = lammps_pot_dihed_nvar(tc);
            d->pot_dihed = lammps_pot_read(tp,tc,1,0,0,nv,d->n_dihed_types,f);
            break;
          case LAMMPS_POT_PAIR:
            tc = lammps_pot_pair_id(v[3]);
            d->pot_pair = lammps_pot_read(tp,tc,1,0,0,2,d->n_atom_types,f);
            break;
          default:
            msg_error_f("\"%s\" potential is not supported",1,v[0]);
            break;
          }
        continue;
        }
      }
    if (n>=3) {
      /* number of atom types */
      if (str_compare(v[1],"atom") && 
          str_compare(v[2],"types")) {
        if (sscanf(v[0],"%u",&(d->n_atom_types))!=1)
          msg_error("invalid specification of number of atom types",1);
        continue;
        }
      /* number of bond types */
      else if (str_compare(v[1],"bond") && 
          str_compare(v[2],"types")) {
        if (sscanf(v[0],"%u",&(d->n_bond_types))!=1)
          msg_error("invalid specification of number of bond types",1);
        continue;
        }
      /* number of angle types */
      else if (str_compare(v[1],"angle") && 
          str_compare(v[2],"types")) {
        if (sscanf(v[0],"%u",&(d->n_angle_types))!=1)
          msg_error("invalid specification of number of angle types",1);
        continue;
        }
      /* number of dihedral types */
      else if (str_compare(v[1],"dihedral") && 
          str_compare(v[2],"types")) {
        if (sscanf(v[0],"%u",&(d->n_dihed_types))!=1)
          msg_error("invalid specification of number of dihedral types",1);
        continue;
        }
      /* number of improper types */
      else if (str_compare(v[1],"improper") && 
          str_compare(v[2],"types")) {
        if (sscanf(v[0],"%u",&(d->n_improp_types))!=1)
          msg_error("invalid specification of number of improper types",1);
        continue;
        }
      /* atoms */
      else if (str_compare(v[0],"atoms")) {
        if (!str_compare(v[1],"#"))
          msg_error("atom style not specified",1);
        d->atom_style = lammps_atom_type_id(v[2]);
        if (!d->atom_style)
          msg_error_f("unknown atom style \"%s\"",1,v[2]);
        nv = lammps_atom_type_nvar(d->atom_style);
        d->atom = lammps_atom_read(d->atom_style,nv,d->n_atoms,f);
        continue;
        }
      }
    else if (n>=2) {
      /* number of atoms */
      if (str_compare(v[1],"atoms")) {
        if (sscanf(v[0],"%u",&(d->n_atoms))!=1)
          msg_error("invalid specification of number of atoms",1);
        continue;
        }
      /* number of bonds */
      else if (str_compare(v[1],"bonds")) {
        if (sscanf(v[0],"%u",&(d->n_bonds))!=1)
          msg_error("invalid specification of number of bonds",1);
        continue;
        }
      /* number of angles */
      else if (str_compare(v[1],"angles")) {
        if (sscanf(v[0],"%u",&(d->n_angles))!=1)
          msg_error("invalid specification of number of angles",1);
        continue;
        }
      /* number of dihedrals */
      else if (str_compare(v[1],"dihedrals")) {
        if (sscanf(v[0],"%u",&(d->n_dihedrals))!=1)
          msg_error("invalid specification of number of dihedrals",1);
        continue;
        }
      /* number of impropers */
      else if (str_compare(v[1],"impropers")) {
        if (sscanf(v[0],"%u",&(d->n_impropers))!=1)
          msg_error("invalid specification of number of impropers",1);
        continue;
        }
      /* number of ellipsoids */
      else if (str_compare(v[1],"ellipsoids")) {
        if (sscanf(v[0],"%u",&(d->n_ellipsoids))!=1)
          msg_error("invalid specification of number of ellipsoids",1);
        continue;
        }
      /* number of lines */
      else if (str_compare(v[1],"lines")) {
        if (sscanf(v[0],"%u",&(d->n_lines))!=1)
          msg_error("invalid specification of number of lines",1);
        continue;
        }
      /* number of triangles */
      else if (str_compare(v[1],"triangles")) {
        if (sscanf(v[0],"%u",&(d->n_triangles))!=1)
          msg_error("invalid specification of number of triangles",1);
        continue;
        }
      /* number of bodies */
      else if (str_compare(v[1],"bodies")) {
        if (sscanf(v[0],"%u",&(d->n_bodies))!=1)
          msg_error("invalid specification of number of bodies",1);
        continue;
        }
      }
    else if (n>=1) {
      /* atomic data */
      if (str_compare(v[0],"atoms")) {
        d->atom_style = LAMMPS_ATOM_CHRG;
        nv = lammps_atom_type_nvar(d->atom_style);
        d->atom = lammps_atom_read(d->atom_style,nv,d->n_atoms,f);
        continue;
        }
      /* particle types */
      if (str_compare(v[0],"masses")) {
        d->type = lammps_type_read(d->n_atom_types,f);
        continue;
        }
      /* atomic velocities */
      else if (str_compare(v[0],"velocities")) {
        nv = lammps_atom_vel_nvar(d->atom_style);
        lammps_atom_read_vel(d->atom,d->atom_style,nv,d->n_atoms,f);
        continue;
        }
      /* interatomic bonds */
      else if (str_compare(v[0],"bonds")) {
        d->bond = lammps_pot_read(LAMMPS_POT_BOND,0,1,1,2,0,d->n_bonds,f);
        continue;
        }
      /* valence angles */
      else if (str_compare(v[0],"angles")) {
        d->angle = lammps_pot_read(LAMMPS_POT_ANGL,0,1,1,3,0,d->n_angles,f);
        continue;
        }
      /* dihedral angles */
      else if (str_compare(v[0],"dihedrals")) {
        d->dihed = lammps_pot_read(LAMMPS_POT_DIHE,0,1,1,4,0,d->n_dihedrals,f);
        continue;
        }
      /* improper angles */
      else if (str_compare(v[0],"impropers")) {
        d->impr = lammps_pot_read(LAMMPS_POT_IMPR,0,1,1,4,0,d->n_impropers,f);
        continue;
        }
      }
    vec_sfree(v,n);
    msg_error_f("unknown data line: %s\n",1,line);
    }
  /* cell vectors */
  lammps_box_set_vectors(d->box);
  /* close the file */
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
