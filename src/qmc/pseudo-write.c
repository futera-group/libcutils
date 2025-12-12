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
#include <cmn/print.h>
#include <mol/atom.h>
#include "qmc/pseudo.h"

/* -------------------------------------------------------------------------- */

/* write pseudopotential data into the file - atom section

   p - pointer to pseudopotential data struct
   f - pointer to open output file */
void pseudo_write_atom(struct pseudo *p, FILE *f) {
  fprintf(f,"&ATOM\n");
  fprintf(f," Z  =%5d\n",p->atom_num);
  fprintf(f," ZV =%5d\n",p->val_num);
  fprintf(f," XC =%5d%15.6f\n",p->xc_num,p->xc_slater);
  switch (p->pp_type) {
    case PSEUDO_NC_NUM: fprintf(f," TYPE = NORMCONSERVING NUMERIC\n"); break;
    default: fprintf(f," TYPE = UNKNOWN\n"); break;
    }
  fprintf(f,"&END\n");
  }

/* write pseudopotential data into the file - info section

   p - pointer to pseudopotential data struct
   f - pointer to open output file */
void pseudo_write_info(struct pseudo *p, FILE *f) {
  unsigned i;
  fprintf(f,"&INFO\n");
  fprintf(f,"    ");
  print_fHline(f,'=',60);
  fprintf(f,"%5s%26s%29s%4s\n","|","Pseudopotential Report",p->time_stamp,"|");
  fprintf(f,"    ");
  print_fHline(f,'-',60);
  fprintf(f,"%5s  %-32s:   %-3s%18s\n","|","Atomic Symbol",
    atom_name(p->atom_num),"|");
  fprintf(f,"%5s  %-32s:   %-3d%18s\n","|","Atomic Number",p->atom_num,"|");
  fprintf(f,"%5s  %-32s:   %-3d%18s\n","|","Number of core states",
    p->state_core,"|");
  fprintf(f,"%5s  %-32s:   %-3d%18s\n","|","Number of valence states",
    p->state_val,"|");
  fprintf(f,"%5s  %-32s:   %21s\n","|","Exchange-Correlation Functional","|");
  fprintf(f,"%5s     %-16s: %7.4f%29s\n","|","Slater exchange",
    p->xc_slater,"|");
  fprintf(f,"%5s     %-16s: %-35s|\n","|","LDA correlation",
    pseudo_func_name(p->xc_lda_c));
  fprintf(f,"%5s     %-16s: %-35s|\n","|","Exchange GC",
    pseudo_func_name(p->xc_gc_e));
  fprintf(f,"%5s     %-16s: %-35s|\n","|","Correlation GC",
    pseudo_func_name(p->xc_gc_c));
  fprintf(f,"%5s  %-23s:%4s%4s%12s%13s\n","|","Electron Configuration",
    "N","L","Occupation","|");
  for (i=0; i<(p->state_core+p->state_val); i++)
    fprintf(f,"%5s%30d%4s%10.4f%15s\n","|",
      p->el_conf_n[i],pseudo_shell_name(p->el_conf_l[i]),
      p->el_conf_occ[i],"|");
  fprintf(f,"%5s  %s%13.6f%17s\n","|","Full Potential Total Energy",
    p->en_pot,"|");
  fprintf(f,"%5s  %s%22s\n","|","Trouiller-Martins normconserving PP","|");
  fprintf(f,"%5s  %4s%5s%10s%13s%25s\n","|","n","l","rc","energy","|");
  for (i=0; i<p->mt_num; i++)
    fprintf(f,"%5s  %4d%5s%10.4f%13.5f%25s\n","|",p->mt_pp_n[i],
      pseudo_shell_name(p->mt_pp_l[i]),p->mt_pp_r[i],p->mt_pp_e[i],"|");
  fprintf(f,"%5s  %s%6d%28s\n","|","Number of Mesh Points :",
    p->mesh_num,"|");
  fprintf(f,"%5s  %s%12.6f%22s\n","|","Pseudoatom Total Energy",
    p->en_tot,"|");
  fprintf(f,"    ");
  print_fHline(f,'=',60);
  fprintf(f,"&END\n");
  }

/* write pseudopotential data into the file - potential section

   p - pointer to pseudopotential data struct
   f - pointer to open output file */
void pseudo_write_potential(struct pseudo *p, FILE *f) {
  unsigned i,j;
  fprintf(f,"&POTENTIAL\n");
  fprintf(f,"%6d\n",p->mesh_num);
  for (i=0; i<p->mesh_num; i++) {
    fprintf(f,"%15.7e",p->dat_pot[i][0]);
    if (p->mt_num)
      fprintf(f,"     ");
    for (j=1; j<(p->mt_num+1); j++)
      fprintf(f,"%15.7e",p->dat_pot[i][j]);
    fprintf(f,"\n");
    }
  fprintf(f,"&END\n");
  }

/* write pseudopotential data into the file - wavefunction section

   p - pointer to pseudopotential data struct
   f - pointer to open output file */
void pseudo_write_wavefunction(struct pseudo *p, FILE *f) {
  unsigned i,j;
  fprintf(f,"&WAVEFUNCTION\n");
  fprintf(f,"%6d CHANNELS=1\n",p->mesh_num);
  for (i=0; i<p->mesh_num; i++) {
    fprintf(f,"%15.7e",p->dat_wfce[i][0]);
    if (p->mt_num)
      fprintf(f,"     ");
    for (j=1; j<(p->mt_num+1); j++)
      fprintf(f,"%15.7e",p->dat_wfce[i][j]);
    fprintf(f,"\n");
    }
  fprintf(f,"&END\n");
  }

/* write pseudopotential data into the file - atom density section

   p - pointer to pseudopotential data struct
   f - pointer to open output file */
void pseudo_write_density(struct pseudo *p, FILE *f) {
  unsigned i,j;
  if (!p->dat_den)
    return;
  fprintf(f,"&ATDENS\n");
  fprintf(f,"%6d\n",p->mesh_num);
  for (i=0; i<p->mesh_num; i++) {
    fprintf(f,"%15.7e",p->dat_wfce[i][0]);
    if (p->mt_num)
      fprintf(f,"     ");
    for (j=1; j<2; j++)
      fprintf(f,"%15.7e",p->dat_wfce[i][j]);
    fprintf(f,"\n");
    }
  fprintf(f,"&END\n");
  }

/* write pseudopotential data into the file

   p - pointer to pseudopotential data struct
   f - name of the output file */
void pseudo_write(struct pseudo *p, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  pseudo_write_atom(p,file);
  pseudo_write_info(p,file);
  pseudo_write_potential(p,file);
  pseudo_write_wavefunction(p,file);
  pseudo_write_density(p,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
