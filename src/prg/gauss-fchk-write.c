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
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* write gaussian data to formatted checkpoint file
 
   g - pointer to gaussian data file
   f - name of output file */
void gauss_fchk_write(struct gauss_dat *g, char *f) {
  unsigned i,j,k,id,n_shells;
  FILE *file = stdout;
  n_shells = basis_shell_num(g->bs);
  if (f && f[0])
    file = file_open(f,"w");
  fprintf(file,"%s\n",g->job_title);
  fprintf(file,"%-10s %-60s %s\n",gauss_job_type_name(g->job_type_id),
    g->job_mthd,g->job_basis);
  fprintf(file,"%-43s%-4s%14d\n","Number of atoms","I",g->n_atoms);
  fprintf(file,"%-43s%-4s%14d\n","Charge","I",g->charge);
  fprintf(file,"%-43s%-4s%14d\n","Multiplicity","I",g->multiplicity);
  fprintf(file,"%-43s%-4s%14d\n","Number of electrons","I",g->n_electrons);
  fprintf(file,"%-43s%-4s%14d\n","Number of alpha electrons","I",
    g->n_alpha_electrons);
  fprintf(file,"%-43s%-4s%14d\n","Number of beta electrons","I",
    g->n_beta_electrons);
  fprintf(file,"%-43s%-4s%14d\n","Number of basis functions","I",
    g->bs->n_bfce);
  fprintf(file,"%-43s%-4s%14d\n","Number of independent functions","I",
    g->bs->n_ibfce);
  fprintf(file,"%-43s%-4sN=%12d\n","Atomic numbers","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].num);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if (g->n_atoms%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Nuclear charges","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].charge[GAUSS_CH_NC]);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if (g->n_atoms%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Current cartesian coordinates","R",
    3*g->n_atoms);
  for (i=0; i<g->n_atoms; i++)
    for (j=0; j<3; j++) {
      fprintf(file,"%16.8e",g->atom[i].coord[j]);
      if (((3*i+j+1)%5)==0)
        fprintf(file,"\n");
      }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Int Atom Types","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].type_id[0]);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Integer atomic weights","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].weight);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Atom fragment info","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].frag_info);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Atom residue num","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].res_info);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Nuclear spins","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].nuc_spin);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Nuclear QMom","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].nuc_qmom);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if (g->n_atoms%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Nuclear GFac","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].nuc_gfac);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if (g->n_atoms%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4s%14d\n","Number of contracted shells","I",
    g->bs->n_cont_shells);
  fprintf(file,"%-43s%-4s%14d\n","Number of primitive shells","I",
    g->bs->n_prim_shells);
  fprintf(file,"%-43s%-4s%14d\n","Pure/Cartesian d shells","I",
    g->bs->n_pure_d);
  fprintf(file,"%-43s%-4s%14d\n","Pure/Cartesian f shells","I",
    g->bs->n_pure_f);
  fprintf(file,"%-43s%-4s%14d\n","Highest angular momentum","I",
    g->bs->max_ang_mom);
  fprintf(file,"%-43s%-4s%14d\n","Largest degree of contraction","I",
    g->bs->max_bfce_cont);
  fprintf(file,"%-43s%-4sN=%12d\n","Shell types","I",n_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++) {
      fprintf(file,"%12d",gauss_bs_shell_id(g->bs->center[i].shell[j].type));
      if (((++id)%6)==0)
        fprintf(file,"\n");
      }
  if (n_shells%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Number of primitives per shell","I",
    n_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++) {
      fprintf(file,"%12d",g->bs->center[i].shell[j].n_prim);
      if (((++id)%6)==0)
        fprintf(file,"\n");
      }
  if (n_shells%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Shell to atom map","I",n_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++) {
      fprintf(file,"%12d",i+1);
      if (((++id)%6)==0)
        fprintf(file,"\n");
      }
  if (n_shells%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Primitive exponents","R",
    g->bs->n_prim_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++)
      for (k=0; k<g->bs->center[i].shell[j].n_prim; k++) {
        fprintf(file,"%16.8e",g->bs->center[i].shell[j].exp[k]);
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
  if (g->bs->n_prim_shells%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Contraction coefficients","R",
    g->bs->n_prim_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++)
      for (k=0; k<g->bs->center[i].shell[j].n_prim; k++) {
        fprintf(file,"%16.8e",g->bs->center[i].shell[j].cf1[k]);
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
  if (g->bs->n_prim_shells%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","P(S=P) Contraction coefficients","R",
    g->bs->n_prim_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++)
      for (k=0; k<g->bs->center[i].shell[j].n_prim; k++) {
        fprintf(file,"%16.8e",(g->bs->center[i].shell[j].type==BASIS_SHELL_SP ?
          g->bs->center[i].shell[j].cf2[k] : 0.0));
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
  if (g->bs->n_prim_shells%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Coordinates of each shell","R",3*n_shells);
  id = 0;
  for (i=0; i<g->bs->n_centers; i++)
    for (j=0; j<g->bs->center[i].n_shells; j++)
      for (k=0; k<3; k++) {
        fprintf(file,"%16.8e",g->bs->center[i].coord[k]);
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
  if ((3*n_shells)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","MM charges","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].charge[GAUSS_CH_MM]);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Real atomic weights","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].mass);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Nuclear ZEff","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].charge[GAUSS_CH_NEFF]);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Constraint Structure","R",
    3*g->n_atoms);
  for (i=0; i<g->n_atoms; i++)
    for (j=0; j<3; j++) {
      fprintf(file,"%16.8e",g->atom[i].coord_fix[j]);
      if (((3*i+j+1)%5)==0)
        fprintf(file,"\n");
      }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  if (g->mo_a) {
    fprintf(file,"%-43s%-4sN=%12d\n","Alpha Orbital Energies","R",
      g->bs->n_ibfce);
    for (i=0; i<g->bs->n_ibfce; i++) {
      fprintf(file,"%16.8e",g->mo_a[i].energy);
      if (((i+1)%5)==0)
        fprintf(file,"\n");
      }
    if (g->bs->n_ibfce%5)
      fprintf(file,"\n");
    }
  if (g->mo_b) {
    fprintf(file,"%-43s%-4sN=%12d\n","Beta Orbital Energies","R",
      g->bs->n_ibfce);
    for (i=0; i<g->bs->n_ibfce; i++) {
      fprintf(file,"%16.8e",g->mo_b[i].energy);
      if (((i+1)%5)==0)
        fprintf(file,"\n");
      }
    if (g->bs->n_ibfce%5)
      fprintf(file,"\n");
    }
  if (g->mo_a) {
    fprintf(file,"%-43s%-4sN=%12d\n","Alpha MO coefficients","R",
      g->bs->n_ibfce*g->bs->n_bfce);
    id = 0;
    for (i=0; i<g->bs->n_ibfce; i++)
      for (j=0; j<g->bs->n_bfce; j++) {
        fprintf(file,"%16.8e",g->mo_a[i].coeff[j]);
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
    if ((g->bs->n_ibfce*g->bs->n_bfce)%5)
      fprintf(file,"\n");
    }
  if (g->mo_b) {
    fprintf(file,"%-43s%-4sN=%12d\n","Beta MO coefficients","R",
      g->bs->n_ibfce*g->bs->n_bfce);
    id = 0;
    for (i=0; i<g->bs->n_ibfce; i++)
      for (j=0; j<g->bs->n_bfce; j++) {
        fprintf(file,"%16.8e",g->mo_b[i].coeff[j]);
        if (((++id)%5)==0)
          fprintf(file,"\n");
        }
    if ((g->bs->n_ibfce*g->bs->n_bfce)%5)
      fprintf(file,"\n");
    }
  fprintf(file,"%-43s%-4sN=%12d\n","Cartesian Gradient","R",
    3*g->n_atoms);
  for (i=0; i<g->n_atoms; i++)
    for (j=0; j<3; j++) {
      fprintf(file,"%16.8e",g->atom[i].grad[j]);
      if (((3*i+j+1)%5)==0)
        fprintf(file,"\n");
      }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Mulliken Charges","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].charge[GAUSS_CH_MULL]);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Int Atom Modified Types","I",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%12d",g->atom[i].type_id[1]);
    if (((i+1)%6)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%6)
    fprintf(file,"\n");
  fprintf(file,"%-43s%-4sN=%12d\n","Atom Modified MM Charges","R",g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    fprintf(file,"%16.8e",g->atom[i].charge[GAUSS_CH_MM_MOD]);
    if (((i+1)%5)==0)
      fprintf(file,"\n");
    }
  if ((3*g->n_atoms)%5)
    fprintf(file,"\n");
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
