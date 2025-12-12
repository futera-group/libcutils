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
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* print out molecular orbital coefficients to file

   g - pointer to gaussian data struct
   a - spin type
   p - occupied only
   f - pointer to open output stream */
void gauss_mol_mo_fprint(struct gauss_dat *g, short a, short p, FILE *f) {
  unsigned ne = (a==GAUSS_SPIN_A ? g->n_alpha_electrons : g->n_beta_electrons);
  unsigned i,ic,ib,id,im,ip,is,nb,nt,nm,bw,fn;
  char atom[80],tp = (a==GAUSS_SPIN_A ? 'A' : 'B');
  struct gauss_mo *m = (a==GAUSS_SPIN_A ? g->mo_a : g->mo_b);
  struct basis_center *c;
  struct basis_shell *s;
  if (!g->bs || !m)
    return;
  /* initialization */
  bw = 5;
  nt = (p ? ne : g->bs->n_ibfce);
  nb = (nt%bw ? nt/bw+1 : nt/bw);
  id = 0;
  /* print out */
  for (ib=0; ib<nb; ib++) {
    nm = ((ib==nb-1) ? (nt%bw ? nt%bw : bw) : bw);
    /* header */
    if (ib)
      fprintf(f,"\n");
    fprintf(f,"              ");
    for (im=0; im<nm; im++)
      fprintf(f,"   %c%05d-%s",tp,id+im+1,((id+im)<ne ? "Occ" : "Vir"));
    fprintf(f,"\n");
    fprintf(f,"              ");
    for (i=0; i<(13*nm); i++)
      fprintf(f,"%c",'-');
    fprintf(f,"\n");
    /* coefficients */
    for (ic=0; ic<g->bs->n_centers; ic++) {
      c = g->bs->center+ic;
      for (is=0; is<c->n_shells; is++) {
        s = c->shell+is;
        for (ip=0; ip<s->n_bfce; ip++) {
          fn = s->bf+ip;
          sprintf(atom,"%s%d",atom_name(c->type),ic+1);
          fprintf(f,"%5d %-4.4s %-3.3s",
            fn+1,atom,basis_shell_name_fce(s->type,ip));
          for (im=0; im<nm; im++)
            fprintf(f,"%13.7f",m[id+im].coeff[fn]);
          fprintf(f,"\n");
          }
        }
      }
    id += nm;
    }
  }

/* print out molecular orbital coefficients

   d - pointer to gaussian data struct
   s - spin type
   c - occupied only
   f - name of the output file */
void gauss_mol_mo_print(struct gauss_dat *d, short s, short c, char *f) {
  FILE *file = stdout;
  /* open file */
  if (f && f[0])
    file = file_open(f,"w");
  /* print coordinates */
  gauss_mol_mo_fprint(d,s,c,file);
  /* close file */
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
