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
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write molecular data to file in apc format
 
   p     - amber prep data structure
   first - first-residuum database indicator
   last  - last-residuum database indicator
   f     - open file stream for output */
void mol_apc_fwrite(struct apc_mol *p, short first, short last, FILE *f) {
  short sgn[3][3] = {{1,1,-1},{1,-1,1},{-1,1,1}};
  char tree_mark[10],s[4][80];
  unsigned i,j;
  /* header */
  if (first)
    fprintf(f,"    0    0    2\n\n");
  fprintf(f,"%s\n",p->title);
  fprintf(f,"%s\n",p->file);
  fprintf(f,"%-3s   XYZ  0\n",p->resname);
  fprintf(f,"CHANGE     OMIT DU   BEG\n");
  fprintf(f,"%8.4f\n",0.0);
  for (i=0; i<3; i++)
    fprintf(f,"%4d  DUMM  DU    M%15.3f%12.3f%12.3f%15.3f\n",
      i+1,999.0*sgn[i][0],999.0*sgn[i][1],999.0*sgn[i][2],0.0);
  /* atom specification */
  for (i=0; i<p->n_atoms; i++) {
    mol_apc_tree_mark(p->atom[i].tree,tree_mark);
    fprintf(f,"%4d  %-4s  %-6s%-4s%15.6f%12.6f%12.6f%12.6f\n",
      p->atom[i].id,p->atom[i].name,p->atom[i].type,tree_mark,
      p->atom[i].coord[0],p->atom[i].coord[1],p->atom[i].coord[2],
      p->atom[i].charge);
    }
  fprintf(f,"\n");
  /* loop specification */
  if (p->n_loops) {
    fprintf(f,"\nLOOP\n");
    for (i=0; i<p->n_loops; i++) {
      for (j=0; j<2; j++)
        str_trim_copy(p->atom[p->loop[i][j]].name,s[j]);
      fprintf(f,"%5s%5s\n",s[0],s[1]);
      }
    fprintf(f,"\n");
    }
  /* improper specification */
  if (p->n_imprs) {
    fprintf(f,"\nIMPROPER\n");
    for (i=0; i<p->n_imprs; i++) {
      for (j=0; j<4; j++)
        str_trim_copy(p->atom[p->impr[i][j]].name,s[j]);
      fprintf(f,"%5s%5s%5s%5s\n",s[0],s[1],s[2],s[3]);
      }
    fprintf(f,"\n");
    }
  fprintf(f,"DONE\n");
  if (last)
    fprintf(f,"STOP\n");
  }

/* write molecular data to file in apc format
 
   p     - amber prep data structure
   first - first-residuum database indicator
   last  - last-residuum database indicator
   f     - name of amber prep file */
void mol_apc_write(struct apc_mol *p, short first, short last, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  mol_apc_fwrite(p,first,last,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
