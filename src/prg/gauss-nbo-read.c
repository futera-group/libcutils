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

#include <math.h>
#include <stdio.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/units.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* read data from NBO output files (basis set file + data file)
 
   d  - pointer to gaussian data struct
   fb - NBO basis set file
   fd - NBO MO data file */
void gauss_nbo_read(struct gauss_dat *g, char *fb, char *fd) {
  char *line;
  unsigned i,j,k,n,nr,nv;
  FILE *file;
  /* basis set */
  g->bs = basis_new();
  basis_nbo_read(g->bs,fb);
  /* atoms */
  g->n_atoms = g->bs->n_centers;
  g->atom = gauss_atom_new(g->n_atoms);
  for (i=0; i<g->n_atoms; i++) {
    for (j=0; j<3; j++)
      g->atom[i].coord[j] = g->bs->center[i].coord[j]/CONV_B_ANG;
    g->atom[i].num = g->bs->center[i].type;
    g->atom[i].charge[GAUSS_CH_NC] = (double)(g->atom[i].num);
    g->atom[i].mass = atom_mass(g->atom[i].num);
    g->atom[i].weight = (unsigned)round(g->atom[i].mass);
    }
  /* MO coefficients */
  file = file_open(fd,"r");
  for (i=0; i<3; i++) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of \"%s\" file while reading header",1,fd);
    line = str_free(line);
    }
  g->mo_a = gauss_mo_new(g->bs->n_ibfce);
  nr = (g->bs->n_ibfce%5 ? g->bs->n_ibfce/5+1 : g->bs->n_ibfce);
  for (i=0; i<g->bs->n_ibfce; i++) {
    n = 0;
    g->mo_a[i].coeff = vec_falloc(g->bs->n_ibfce);
    for (j=0; j<nr; j++) {
      line = str_read_line_new(file);
      if (!line)
        msg_error("invalid format of MO coefficient data array",1);
      nv = ((j+1)<nr ? 5 : (g->bs->n_ibfce%5 ? g->bs->n_ibfce%5 : 5));
      for (k=0; k<nv; k++)
        if (sscanf(line+1+15*k,"%lf",&(g->mo_a[i].coeff[n++]))!=1)
          msg_error("cannot read MO coefficient data array",1);
      line = str_free(line);
      }
    }
  file_close(file);
  /* reset basis coordinates */
  gauss_bs_set_coord(g);
  }

/* -------------------------------------------------------------------------- */
