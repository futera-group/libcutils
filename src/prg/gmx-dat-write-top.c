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
#include <cmn/list.h>
#include <cmn/string.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Write force-field parameters to open gromacs topology file 
 
   p - force field data
   f - open file stream */
void gmx_dat_fwrite_top_ff(struct gmx_ff *p, FILE *f) {
  unsigned i,j;
  struct gmx_atom *a;
  struct gmx_bond *b;
  /* default setting */
  fprintf(f,"; default force-field settings\n");
  fprintf(f,"[ defaults ]\n");
  fprintf(f,"; %5s%15s%15s%15s%15s\n","nb-func","comb-rule",
    "gen-pairs","scale-lj","scale-qq");
  fprintf(f,"%6d%14d%15s%19.6f%15.6f\n",p->nb_func,p->comb_rule,
    p->gen_pairs ? "yes" : "no",p->scale_lj,p->scale_qq);
  fprintf(f,"\n");
  /* atom types */
  if (p->dat->n_atom_types) {
    fprintf(f,"; force-field atomic types\n");
    fprintf(f,"[ atomtypes ]\n");
    fprintf(f,"; %4s%6s%12s%17s%9s%10s%18s\n","type","num",
      "mass","charge","particle","sigma","epsilon");
    for (i=0; i<p->dat->n_atom_types; i++) {
      a = p->dat->type + i;
      fprintf(f,"  %-5s%5d%15.6f%15.6f%5s%16.10f%15.6f\n",
        a->type,a->num,a->mass,a->charge,gmx_ff_particle_type_name(a->particle),
        a->sigma,a->epsilon);
      }
    fprintf(f,"\n");
    }
  /* pair types */
  if (p->dat->n_14_pairs) {
    fprintf(f,"; force-field 1-4 pair types\n");
    fprintf(f,"[ pairtypes ]\n");
    fprintf(f,"; %2s%5s%9s%14s%19s\n","i","j","func","sigma","epsilon");
    for (i=0; i<p->dat->n_14_pairs; i++) {
      b = p->dat->pair + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%18.10f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* non-bonded parameters */
  if (p->dat->n_nb_pairs) {
    fprintf(f,"; force-field other non-bonded pair types\n");
    fprintf(f,"[ nonbond_params ]\n");
    fprintf(f,"; %2s%5s%9s%14s%19s\n","i","j","func","sigma","epsilon");
    for (i=0; i<p->dat->n_nb_pairs; i++) {
      b = p->dat->nbpr + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%18.10f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* bond types */
  if (p->dat->n_bonds) {
    fprintf(f,"; force-field bond types\n");
    fprintf(f,"[ bondtypes ]\n");
    fprintf(f,"; %2s%5s%9s%12s%12s\n","i","j","func","b0","kb");
    for (i=0; i<p->dat->n_bonds; i++) {
      b = p->dat->bond + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.4f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* constrain types */
  if (p->dat->n_constraints) {
    fprintf(f,"; force-field constraint types\n");
    fprintf(f,"[ constrainttypes ]\n");
    fprintf(f,"; %2s%6s%10s%11s\n","i","j","func","c0");
    for (i=0; i<p->dat->n_constraints; i++) {
      b = p->dat->cnst + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-6s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.6f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* angle types */
  if (p->dat->n_angles) {
    fprintf(f,"; force-field angle types\n");
    fprintf(f,"[ angletypes ]\n");
    fprintf(f,"; %2s%5s%5s%9s%12s%15s%15s%15s\n","i","j","k",
      "func","th0","cth","ub0","cub");
    for (i=0; i<p->dat->n_angles; i++) {
      b = p->dat->angl + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.4f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* proper dihedral types */
  if (p->dat->n_dihedrals) {
    fprintf(f,"; force-field proper dihedral types\n");
    fprintf(f,"[ dihedraltypes ]\n");
    fprintf(f,"; %2s%5s%5s%5s%9s%12s%15s%15s\n","i","j","k","l",
      "func","phi0","cp","mult");
    for (i=0; i<p->dat->n_dihedrals; i++) {
      b = p->dat->dihe + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.4f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* improper dihedral types */
  if (p->dat->n_impropers) {
    fprintf(f,"; force-field improper dihedral types\n");
    fprintf(f,"[ dihedraltypes ]\n");
    fprintf(f,"; %2s%5s%5s%5s%9s%12s%15s\n","i","j","k","l",
      "func","q0","cq");
    for (i=0; i<p->dat->n_impropers; i++) {
      b = p->dat->impr + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.4f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* implicit solvent parameters */
  if (p->dat->n_igb_parms) {
    fprintf(f,"; force-field implicit solvent parameters\n");
    fprintf(f,"[ implicit_genborn_params ]\n");
    fprintf(f,"; %-10s%8s%14s%15s%16s%15s\n","type",
      "sar","st","pi","gbr","hct");
    for (i=0; i<p->dat->n_igb_parms; i++) {
      b = p->dat->igbp + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      for (j=0; j<b->n_parms; j++)
        fprintf(f,"%15.6f",b->parm[j]);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  /* cmap types */
  if (p->dat->n_cmap) {
    fprintf(f,"; force-field cmap types\n");
    fprintf(f,"[ cmaptypes ]\n");
    fprintf(f,"; %2s%5s%5s%5s%5s\n","i","j","k","l","n");
    for (i=0; i<p->dat->n_cmap; i++) {
      b = p->dat->cmap + i;
      fprintf(f,"  ");
      for (j=0; j<b->n_atoms; j++)
        fprintf(f,"%-5s",b->name[j]);
      fprintf(f,"%5d",b->type);
      if (b->n_parms)
        fprintf(f,"%5d%5d %c\n",(int)b->parm[0],(int)b->parm[1],'\\');
      for (j=2; j<b->n_parms; j++) {
        fprintf(f,"%16.10f",b->parm[j]);
        if ((j+1)<b->n_parms && (j-1)%6==0)
          fprintf(f," %c\n",'\\');
        }
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  }

/* Write structure to open gromacs topology file

   d    - gromacs data structure
   full - write explicitly all data rather than using inclusion statements
   chrg - atomic charges to rewrite saved data (not used if NULL)
   f    - open file stream */
void gmx_dat_fwrite_top(struct gmx_dat *d, short full, double *chrg, FILE *f) {
  unsigned i,j,id;
  char *name;
  FILE *f_itp;
  struct list *t;
  struct ldata *p;
  struct gmx_mol *m;
  fprintf(f,"\n");
  /* force-field parameters */
  if (!full && d->ff->file && d->ff->file[0]) {
    fprintf(f,"; Force-field parameters\n");
    fprintf(f,"#include \"%s\"\n",d->ff->file);
    fprintf(f,"\n");
    }
  else
    gmx_dat_fwrite_top_ff(d->ff,f);
  /* molecule types */
  t = list_alloc();
  for (i=0; i<d->n_mols; i++)
    if (d->mol[i].file && d->mol[i].file[0])
      if (!list_find_data_t(t,d->mol+i,gmx_mol_cmp_file_wrp_search))
        list_add_end_p(t,d->mol+i);
  /* include ITP if available */
  if ((!full || !chrg) && t->num) {
    fprintf(f,"; Molecular topologies\n");
    list_sort(t,gmx_mol_cmp_file_wrp_sort);
    for (p=t->first; p; p=p->l_next) {
      m = (struct gmx_mol*)p->l_data;
      fprintf(f,"#include \"%s\"\n",m->file);
      if (m->gen_itp) {
        f_itp = file_open(m->file,"w");
        gmx_mol_fwrite_itp(m,NULL,f_itp);
        file_close(f_itp);
        }
      }
    for (i=0; i<d->n_mols; i++)
      if (!d->mol[i].file || !d->mol[i].file[0])
        gmx_mol_fwrite_itp(d->mol+i,NULL,f);
    fprintf(f,"\n");
    }
  /* explicit topology of all molecules */
  else if (chrg) {
    for (i=0,id=0; i<d->n_frags; i++)
      for (j=0; j<d->frag[i].n_rep; j++) {
        name = d->frag[i].mol->name;
        d->frag[i].mol->name = str_sprintf("%s_%d_%d",name,i+1,j+1);
        gmx_mol_fwrite_itp(d->frag[i].mol,chrg+id,f);
        str_free(d->frag[i].mol->name);
        d->frag[i].mol->name = name;
        id += d->frag[i].mol->n_atoms;
        }
    }
  /* explicit topology of each unique molecule */
  else {
    for (i=0; i<d->n_mols; i++)
      gmx_mol_fwrite_itp(d->mol+i,NULL,f);
    }
  list_free(t,NULL);
  /* system */
  fprintf(f,"; total system\n");
  fprintf(f,"[ system ]\n");
  fprintf(f,"; title\n");
  fprintf(f,"%s\n",d->title);
  fprintf(f,"\n");
  /* molecules */
  fprintf(f,"; system composition\n");
  fprintf(f,"[ molecules ]\n");
  fprintf(f,"; %-15s%7s\n","name","num");
  if (chrg) {
    for (i=0; i<d->n_frags; i++)
      for (j=0; j<d->frag[i].n_rep; j++) {
        name = str_sprintf("%s_%d_%d",d->frag[i].name,i+1,j+1);
        fprintf(f,"  %-15s %5d\n",name,1);
        str_free(name);
        }
    }
  else {
    for (i=0; i<d->n_frags; i++)
      fprintf(f,"  %-15s %5d\n",d->frag[i].name,d->frag[i].n_rep);
    }
  fprintf(f,"\n");
  }

/* Write structure to gromacs topology file

   d    - gromacs data structure
   full - write explicitly all data rather than using inclusion statements
   chrg - atomic charges to rewrite saved data (not used if NULL)
   name - name of the file */
void gmx_dat_write_top(struct gmx_dat *d, short full, double *chrg,
  char *name) {
  FILE *f = stdout;
  if (name && name[0])
    f = file_open(name,"w");
  gmx_dat_fwrite_top(d,full,chrg,f);
  if (name && name[0])
    file_close(f);
  }

/* -------------------------------------------------------------------------- */
