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
#include <cmn/matrix.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write molecular data to already open cif file 
 
   c - pointer to cif molecular data struct
   f - pointer to open cif file */
void mol_cif_fwrite(struct cif_mol *c, FILE *f) {
  double x[3];
  unsigned i;
  /* header */
  fprintf(f,"data_image0\n");
  fprintf(f,"%-21s%.3f\n","_cell_length_a",
    c->cell->side[0]);
  fprintf(f,"%-21s%.3f\n","_cell_length_b",
    c->cell->side[1]);
  fprintf(f,"%-21s%.3f\n","_cell_length_c",
    c->cell->side[2]);
  fprintf(f,"%-21s%.1f\n","_cell_angle_alpha",
    c->cell->angle[0]);
  fprintf(f,"%-21s%.1f\n","_cell_angle_beta",
    c->cell->angle[1]);
  fprintf(f,"%-21s%.1f\n\n","_cell_angle_gamma",
    c->cell->angle[2]);
  fprintf(f,"%-34s%s\n","_symmetry_space_group_name_H-M",
    c->cell->space_group_name);
  fprintf(f,"%-34s%d\n\n","_symmetry_int_tables_number",
    c->cell->space_group_id);
  fprintf(f,"loop_\n");
  fprintf(f,"  _symmetry_equiv_pos_as_xyz\n");
  fprintf(f,"  'x, y, z'\n\n");
  fprintf(f,"loop_\n");
  fprintf(f,"  _atom_site_label\n");
  fprintf(f,"  _atom_site_occupancy\n");
  fprintf(f,"  _atom_site_fract_x\n");
  fprintf(f,"  _atom_site_fract_y\n");
  fprintf(f,"  _atom_site_fract_z\n");
  fprintf(f,"  _atom_site_B_iso_or_equiv\n");
  fprintf(f,"  _atom_site_thermal_displace_type\n");
  fprintf(f,"  _atom_site_type_symbol\n");
  /* coordinates */
  for (i=0; i<c->n_atoms; i++) {
    if (c->cell->transf) {
      mat_fmult_rvec(x,c->cell->t_r2c_c,c->atom[i].coord,3,3);
      fprintf(f,"  %-5s%10.3f%9.5f%9.5f%9.5f  Biso%8.3f  %s\n",
        c->atom[i].name,1.0,x[0],x[1],x[2],1.0,atom_name(c->atom[i].num));
      }
    else
      fprintf(f,"  %-5s%10.3f%9.5f%9.5f%9.5f  Biso%8.3f  %s\n",
        c->atom[i].name,1.0,c->atom[i].coord[0],c->atom[i].coord[1],
          c->atom[i].coord[2],1.0,atom_name(c->atom[i].num));
    }
  }

/* write molecular data to cif file
 
   c - pointer to cif molecular data struct
   f - name of the file */
void mol_cif_write(struct cif_mol *c, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  mol_cif_fwrite(c,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
