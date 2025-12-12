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
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* write structure coordinates into open file
 
   d - pointer to cpmd data struct
   s - structure specification
   t - file format 
   f - pointer to open file */
void cpmd_dat_write_geom_f(struct cpmd_dat *d, long int s, short t, FILE *f) {
  struct gen_mol *m;
  struct cif_mol *c;
  struct xyz_mol *x;
  unsigned i;
  /* all geometry optimization */
  if (s<-1) {
    if (t==CPMD_FRMT_XYZ) {
      for (i=0; i<d->n_geoms; i++) {
        m = cpmd_dat_conv_mol(d,(long)(i+1));
        if (m && m->n_atoms) {
          x = mol_gen_xyz(m);
          sprintf(x->title,"%20.10f H",d->g_opt[i].energy);
          mol_xyz_fwrite(x,f);
          mol_xyz_free(x);
          mol_gen_free(m);
          }
        }
      }
    }
  /* one specific structure */
  else {
    m = cpmd_dat_conv_mol(d,s);
    if (m && m->n_atoms) {
      switch (t) {
        /* CIF file format */
        case CPMD_FRMT_CIF: 
          c = mol_gen_cif(m);
          mol_cif_fwrite(c,f);
          mol_cif_free(c);
          break;
        /* XYZ file format */
        case CPMD_FRMT_XYZ: 
          x = mol_gen_xyz(m);
          mol_xyz_fwrite(x,f);
          mol_xyz_free(x);
          break;
        }
      mol_gen_free(m);
      }
    }
  }

/* write structure coordinates into file
 
   d - pointer to cpmd data struct
   s - structure specification
   t - file format 
   f - name of the file */
void cpmd_dat_write_geom(struct cpmd_dat *d, long int s, short t, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cpmd_dat_write_geom_f(d,s,t,file);
  if (f && f[0])
    file_close(file);
  }

/* write energy bands into open file
 
   d - pointer to cpmd data struct
   f - pointer to open file */
void cpmd_dat_write_band_f(struct cpmd_dat *d, FILE *f) {
  unsigned i,j;
  /* spin polarized system */
  if (d->lsd_calc) {
    for (i=0; i<d->n_kpoints; i++) {
      fprintf(f,"%5d",i+1);
      for (j=0; j<d->n_alpha_states; j++)
        fprintf(f,"%15.7e",d->state_a[d->n_alpha_states*i+j].energy);
      fprintf(f,"\n");
      }
    fprintf(f,"\n\n");
    for (i=0; i<d->n_kpoints; i++) {
      fprintf(f,"%5d",i+1);
      for (j=0; j<d->n_beta_states; j++)
        fprintf(f,"%15.7e",d->state_b[d->n_beta_states*i+j].energy);
      fprintf(f,"\n");
      }
    }
  /* closed shell system */
  else {
    for (i=0; i<d->n_kpoints; i++) {
      fprintf(f,"%5d",i+1);
      for (j=0; j<d->n_states; j++)
        fprintf(f,"%15.7e",d->state_a[d->n_states*i+j].energy);
      fprintf(f,"\n");
      }
    }
  }

/* write energy bands into file
 
   d - pointer to cpmd data struct
   f - name of the file */
void cpmd_dat_write_band(struct cpmd_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cpmd_dat_write_band_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* print out projection of one set of states to atomic basis set
 
   vs      - array with projection data
   va      - array with atomic data
   vt      - array with types
   n_state - number of states
   n_atom  - number of atoms
   f       - pointer to open cpmd file */
void cpmd_dat_write_proj_one(struct cpmd_state *vs, struct cpmd_atom *va,
  struct cpmd_type *vt, unsigned n_state, unsigned n_atom, FILE *f) {
  unsigned ib,is,ia,it,ir,ix,id,nb,ns;
  nb = (n_state%8 ? n_state/8+1 : n_state);
  for (ib=0; ib<nb; ib++) {
    ns = (ib==nb-1 ? (n_state%8 ? n_state%8 : 8) : 8);
    /* blank line between blocks */
    if (ib)
      fprintf(f,"\n");
    /* header */
    fprintf(f,"%13s","ORBITAL");
    for (is=0; is<ns; is++)
      fprintf(f,(is ? "%8d" : "%7d"),ib*8+is+1);
    fprintf(f,"\n%13s","COMPLETNESS");
    for (is=0; is<ns; is++)
      fprintf(f,(is ? "%8.3f" : "%9.3f"),vs[ib*8+is].compln);
    fprintf(f,"\n%13s","OCCUPATION");
    for (is=0; is<ns; is++)
      fprintf(f,(is ? "%8.3f" : "%9.3f"),vs[ib*8+is].occup);
    fprintf(f,"\n\n");
    /* coefficients */
    for (ia=0; ia<n_atom; ia++) {
      it = va[ia].t_id;
      for (ir=0; ir<vt[it].n_ao_ids; ir++)
        for (ix=0; ix<vt[it].ao_id[ir]; ix++) {
          if (ir)
            fprintf(f,"         %-4s",
              cpmd_bs_sym(vt[it].ao_id[ir],ix));
          else
            fprintf(f,"%3d%4s  %-4s",ia+1,atom_name(vt[it].num),
              cpmd_bs_sym(vt[it].ao_id[ir],ix));
          for (is=0; is<ns; is++)
            fprintf(f,(is ? "%8.3f" : "%9.3f"),vs[ib*8+is].ao_pj[id]);
          fprintf(f,"\n");
          id++;
          }
      }
    id = 0;
    }
  }

/* write state projection into open file
 
   d - pointer to cpmd data struct
   f - pointer to open file */
void cpmd_dat_write_proj_f(struct cpmd_dat *d, FILE *f) {
  /* spin polarized system */
  if (d->lsd_calc) {
    cpmd_dat_write_proj_one(d->state_a, d->atom, d->type,
      d->n_alpha_states,d->n_atoms,f);
    fprintf(f,"\n\n");
    cpmd_dat_write_proj_one(d->state_b, d->atom, d->type,
      d->n_beta_states,d->n_atoms,f);
    }
  /* closed shell system */
  else {
    cpmd_dat_write_proj_one(d->state_a, d->atom, d->type,
      d->n_states,d->n_atoms,f);
    }
  }

/* write state projection into file
 
   d - pointer to cpmd data struct
   f - name of the file */
void cpmd_dat_write_proj(struct cpmd_dat *d, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  cpmd_dat_write_proj_f(d,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
