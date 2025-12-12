
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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Write specified bonding data to open gromacs topology file

   m  - molecular data structure 
   tp - bonding term type
   f  - open file stream */
void gmx_mol_fwrite_itp_bonding(struct gmx_mol *m, short tp, FILE *f) {
  unsigned i,j,k,*n,na;
  short found;
  struct gmx_res *r;
  struct gmx_bond **b;
  na = gmx_bond_key_natoms(tp);
  /* intra-residual terms */
  for (i=0,found=0; i<m->n_resids; i++) {
    b = gmx_top_bond_get(m->res[i].top,tp,&n);
    if (*n) {
      found = 1;
      break;
      }
    }
  if (found) {
    /* header */
    fprintf(f,"; intra-residual %ss\n",gmx_bond_key_desc(tp));
    fprintf(f,"[ %s ]\n",gmx_bond_key_name(tp));
    fprintf(f,";");
    for (i=0; i<na; i++)
      fprintf(f,"%*c",i ? 5 : 4,(char)((int)'i'+i));
    fprintf(f,"%7s%8s%10s%10s%10s\n","func","c0","c1","c2","c3");
    /* residues */  
    for (i=0; i<m->n_resids; i++) {
      r = m->res+i;
      b = gmx_top_bond_get(m->res[i].top,tp,&n);
      if (*n)
        fprintf(f,"; residuum %5d %5s\n",r->id+1,r->name);
      /* bonding */
      for (j=0; j<(*n); j++) {
        for (k=0; k<(*b)[j].n_atoms; k++)
          fprintf(f,"%5d",
            m->res[(*b)[j].ir[k]].atom[(*b)[j].ia[k]].id+1);
        fprintf(f,"%5d",(*b)[j].type);
        for (k=0; k<(*b)[j].n_parms; k++)
          fprintf(f,"%14.3f",(*b)[j].parm[k]);
        fprintf(f,"  ; ");
        for (k=0; k<(*b)[j].n_atoms; k++)
          fprintf(f,"%-5s",
            m->res[(*b)[j].ir[k]].atom[(*b)[j].ia[k]].name);
        fprintf(f,"\n");
        }
      }
    fprintf(f,"\n");
    }
  /* inter-residual terms */
  b = gmx_top_bond_get(m->top,tp,&n);
  if (*n) {
    /* header */
    fprintf(f,"; inter-residual %ss\n",gmx_bond_key_desc(tp));
    fprintf(f,"[ %s ]\n",gmx_bond_key_name(tp));
    fprintf(f,";");
    for (i=0; i<na; i++)
      fprintf(f,"%*c",i ? 5 : 4,(char)((int)'i'+i));
    fprintf(f,"%7s%8s%10s%10s%10s\n","func","c0","c1","c2","c3");
    /* bonding */
    for (i=0; i<(*n); i++) {
      for (j=0; j<(*b)[i].n_atoms; j++)
        fprintf(f,"%5d",
          m->res[(*b)[i].ir[j]].atom[(*b)[i].ia[j]].id+1);
      fprintf(f,"%5d",(*b)[i].type);
      for (j=0; j<(*b)[i].n_parms; j++)
        fprintf(f,"%14.3f",(*b)[i].parm[j]);
      fprintf(f,"  ; ");
      for (j=0; j<(*b)[i].n_atoms; j++)
        fprintf(f,"%3d(%3.3s:%-4.4s)",
          m->res[(*b)[i].ir[j]].id+1,
          m->res[(*b)[i].ir[j]].name,
          m->res[(*b)[i].ir[j]].atom[(*b)[i].ia[j]].name);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  }

/* Write exclusion list to open gromacs topology file

   m  - molecular data structure 
   f  - open file stream */
void gmx_mol_fwrite_itp_exclusions(struct gmx_mol *m, FILE *f) {
  unsigned i,j,k,*n;
  short found;
  struct gmx_res *r;
  struct gmx_bond **b;
  /* intra-residual terms */
  for (i=0,found=0; i<m->n_resids; i++) {
    b = gmx_top_bond_get(m->res[i].top,GMX_BTERM_EXCL,&n);
    if (*n) {
      found = 1;
      break;
      }
    }
  if (found) {
    /* header */
    fprintf(f,"; intra-residual exclusions\n");
    fprintf(f,"[ exclusions ]\n");
    /* residues */  
    for (i=0; i<m->n_resids; i++) {
      r = m->res+i;
      b = gmx_top_bond_get(m->res[i].top,GMX_BTERM_EXCL,&n);
      if (*n)
        fprintf(f,"; residuum %5d %5s\n",r->id+1,r->name);
      /* bonding */
      for (j=0; j<(*n); j++) {
        for (k=0; k<(*b)[j].n_atoms; k++)
          fprintf(f,"%5d",
            m->res[(*b)[j].ir[k]].atom[(*b)[j].ia[k]].id+1);
        fprintf(f,"  ; ");
        for (k=0; k<(*b)[j].n_atoms; k++)
          fprintf(f,"%-5s",
            m->res[(*b)[j].ir[k]].atom[(*b)[j].ia[k]].name);
        fprintf(f,"\n");
        }
      }
    fprintf(f,"\n");
    }
  /* inter-residual terms */
  b = gmx_top_bond_get(m->top,GMX_BTERM_EXCL,&n);
  if (*n) {
    /* header */
    fprintf(f,"; inter-residual exclusions\n");
    fprintf(f,"[ exclusions ]\n");
    for (i=0; i<(*n); i++) {
      for (j=0; j<(*b)[i].n_atoms; j++)
        fprintf(f,"%5d",
          m->res[(*b)[i].ir[j]].atom[(*b)[i].ia[j]].id+1);
      fprintf(f,"  ; ");
      for (j=0; j<(*b)[i].n_atoms; j++)
        fprintf(f,"%3d(%3.3s:%-4.4s)",
          m->res[(*b)[i].ir[j]].id+1,
          m->res[(*b)[i].ir[j]].name,
          m->res[(*b)[i].ir[j]].atom[(*b)[i].ia[j]].name);
      fprintf(f,"\n");
      }
    fprintf(f,"\n");
    }
  }

/* Write molecular type to open gromacs topology file

   m - molecular data structure 
   q - atomic charges (use stored values if NULL)
   f - open file stream */
void gmx_mol_fwrite_itp(struct gmx_mol *m, double *q, FILE *f) {
  unsigned i,j,ia;
  double qatm,qtot;
  struct gmx_res *r;
  struct gmx_atom *a;
  fprintf(f,"; molecular structure\n");
  fprintf(f,"[ moleculetype ]\n");
  fprintf(f,"; %2s%15s\n","name","n-excl");
  fprintf(f,"  %-15s %d\n",m->name,m->n_excl);
  fprintf(f,"\n");
  /* atoms */
  if (m->n_resids) {
    fprintf(f,"; atomic data\n");
    fprintf(f,"[ atoms ]\n");
    fprintf(f,"; %3s%5s%5s%5s%5s%5s%14s%14s\n","id",
      "type","ir","res","atom","cgnr","charge","mass");
    for (i=0,ia=0,qtot=0.0; i<m->n_resids; i++) {
      r = m->res+i;
      if (r->n_atoms)
        fprintf(f,"; residuum %5d %5s\n",r->id+1,r->name);
      for (j=0; j<r->n_atoms; j++,ia++) {
        a = r->atom+j;
        qatm = (q ? q[ia] : a->charge);
        qtot += a->charge;
        fprintf(f,"%5d%5s%5d%5s%5s%5d%15.6f%15.6f ; q %8.3f\n",
          a->id+1,a->type,r->id+1,r->name,a->name,j+1,qatm,a->mass,qtot);
        }
      }
    fprintf(f,"\n");
    }
  /* bonding */
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_BOND,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_CNST,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_STTL,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_PAIR,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_NBPR,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_ANGL,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_DIHE,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_IMPR,f);
  gmx_mol_fwrite_itp_bonding(m,GMX_BTERM_CMAP,f);
  /* exclusions */
  gmx_mol_fwrite_itp_exclusions(m,f);
  }

/* -------------------------------------------------------------------------- */
