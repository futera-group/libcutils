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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* write potential IDs and coefficients to LAMMPS data file 
 
   p      - potential parameters
   wrt_id - write sequential ID
   wrt_tp - write potential ID
   n_pots - number of the terms
   f      - open file stream */
void lammps_pot_write_data(struct lammps_pot *p, short wrt_id, short wrt_tp,
  unsigned n_pots, FILE *f) {
  unsigned i,j;
  for (i=0; i<n_pots; i++) {
    if (wrt_id)
      fprintf(f,"%10d",i+1);
    if (wrt_tp)
      fprintf(f,"%8d",p[i].type+1);
    for (j=0; j<p[i].n_ids; j++)
      fprintf(f,"%10d",p[i].id[j]+1);
    for (j=0; j<p[i].n_coeffs; j++)
      fprintf(f,"%15.6f",p[i].coeff[j]);
    fprintf(f,"\n");
    }
  }

/* write potential specification to LAMMPS data file 
 
   p        - potential parameters
   pot_type - type of the potential
   out_type - output type (IDs,Coeffs,Both)
   n_pots   - number of the terms
   f        - open file stream */
void lammps_pot_write(struct lammps_pot *p, short pot_type, short out_type,
  unsigned n_pots, FILE *f) {
  /* potential type */
  if (p && n_pots) {
    switch (pot_type) {
      case LAMMPS_POT_BOND:
        if (out_type) {
          fprintf(f,"\nBond Coeffs # %s\n\n",
            lammps_pot_bond_sym(p[0].style));
          lammps_pot_write_data(p,1,0,n_pots,f);
          }
        else {
          fprintf(f,"\nBonds\n\n");
          lammps_pot_write_data(p,1,1,n_pots,f);
          }
        break;
      case LAMMPS_POT_ANGL:
        if (out_type) {
          fprintf(f,"\nAngle Coeffs # %s\n\n",
            lammps_pot_angle_sym(p[0].style));
          lammps_pot_write_data(p,1,0,n_pots,f);
          }
        else {
          fprintf(f,"\nAngles\n\n");
          lammps_pot_write_data(p,1,1,n_pots,f);
          }
        break;
      case LAMMPS_POT_DIHE:
        if (out_type) {
          fprintf(f,"\nDihedral Coeffs # %s\n\n",
            lammps_pot_dihed_sym(p[0].style));
          lammps_pot_write_data(p,1,0,n_pots,f);
          }
        else {
          fprintf(f,"\nDihedrals\n\n");
          lammps_pot_write_data(p,1,1,n_pots,f);
          }
        break;
      case LAMMPS_POT_PAIR:
        fprintf(f,"\nPair Coeffs # %s\n\n",
          lammps_pot_pair_sym(p[0].style));
        lammps_pot_write_data(p,1,0,n_pots,f);
        break;
      }
    }
  }

/* -------------------------------------------------------------------------- */
