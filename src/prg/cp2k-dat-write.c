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
#include <mol/atom.h>
#include <mol/molec.h>
#include <qmc/basis.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* write structure coordinates into open file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_dat_write_geom_f(struct cp2k_dat *d, FILE *f) {
  struct gen_mol *m;
  struct xyz_mol *x;
  if (d && d->n_atoms) {
    m = cp2k_dat_conv_mol(d);
    x = mol_gen_xyz(m);
    mol_xyz_fwrite(x,f);
    mol_xyz_free(x);
    mol_gen_free(m);
    }
  }

/* write structure coordinates into file
 
   d - cp2k data struct
   f - name of the file */
void cp2k_dat_write_geom(struct cp2k_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cp2k_dat_write_geom_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* write basis set into open file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_dat_write_basis_f(struct cp2k_dat *d, FILE *f) {
  struct basis *b;
  if (d && d->n_atoms) {
    b = cp2k_dat_conv_basis(d);
    basis_print(b,1,NULL);
    basis_free(b);
    }
  }

/* write basis set into file
 
   d - cp2k data struct
   f - name of the file */
void cp2k_dat_write_basis(struct cp2k_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cp2k_dat_write_basis_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* write overlap matrix into open file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_dat_write_overlap_f(struct cp2k_dat *d, FILE *f) {
  if (d && d->n_sphr_fce && d->ovrl)
    mat_ffprint_Sc(f,d->ovrl,d->n_sphr_fce,5);
  }

/* write overlap matrix into file
 
   d - cp2k data struct
   f - name of the file */
void cp2k_dat_write_overlap(struct cp2k_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cp2k_dat_write_overlap_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* write electronic states into open file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_dat_write_states_f(struct cp2k_dat *d, FILE *f) {
  unsigned i1,i2,ib1,ia,ib,ic,ik,ix,n_blocks,n_cols,n_vals = 4;
  if (d && d->n_sphr_fce && d->kpoint) {
    /* k-points */
    for (ik=0; ik<d->n_kpoints; ik++) {
      /* spin components */
      for (ix=0; ix<d->n_spins; ix++) {
        if (d->n_spins>1)
          printf("K-point #%d, %s states:\n\n",ik+1,ix ? "alpha" : "beta");
        n_blocks = (d->n_states[ix]%n_vals ? 
          d->n_states[ix]/n_vals+1 : d->n_states[ix]/n_vals);
        /* data blocks */
        for (ib1=0; ib1<n_blocks; ib1++) {
          n_cols = ((ib1+1)<n_blocks ? n_vals :
           (d->n_states[ix]%n_vals ? d->n_states[ix]%n_vals : n_vals));
          i1 = 0;
          /* header */
          i2 = ib1*n_vals;
          fprintf(f,"%19s"," ");
          for (ic=0; ic<n_cols; ic++)
            fprintf(f,"%13d",++i2);
          fprintf(f,"\n");
          /* energies */
          i2 = ib1*n_vals;
          fprintf(f,"%23s"," ");
          for (ic=0; ic<n_cols; ic++)
            fprintf(f,"%13.6f",d->kpoint[ik].state[ix][i2++].energy);
          fprintf(f,"\n\n");
          /* occupancies */
          i2 = ib1*n_vals;
          fprintf(f,"%23s"," ");
          for (ic=0; ic<n_cols; ic++)
            fprintf(f,"%13.6f",d->kpoint[ik].state[ix][i2++].occup);
          fprintf(f,"\n\n");
          /* atoms */
          for (ia=0; ia<d->n_atoms; ia++) {
            /* atomic orbitals */
            for (ib=0; ib<d->kind[d->atom[ia].kind].n_sphr_fce; ib++) {
              fprintf(f,"%6d%6d  %-3s%6d",i1+1,ia+1,
                atom_name(d->atom[ia].num),ib+1);
              /* MO coefficients */
              i2 = ib*n_vals;
              for (ic=0; ic<n_cols; ic++)
                fprintf(f,"%13.6f",d->kpoint[ik].state[ix][i2++].ao_pj[i1]);
              fprintf(f,"\n");
              i1++;
              }
            fprintf(f,"\n");
            }
          }
        }
      }
    }
  }

/* write electronic states into file
 
   d - cp2k data struct
   f - name of the file */
void cp2k_dat_write_states(struct cp2k_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cp2k_dat_write_states_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
