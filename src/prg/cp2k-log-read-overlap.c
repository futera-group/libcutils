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

#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read overlap matrix from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_overlap(struct cp2k_dat *d, FILE *f) {
  unsigned i,j,k,l,i1,i2,n_blocks,n_cols,n_vals = 4;
  char *line;
  /* find the overal data block */
  line = str_ffind_b_new(f," OVERLAP MATRIX");
  if (line) {
    line = str_free(line);
    /* initialization */
    if (d->n_sphr_fce && d->n_ibfce!=d->n_sphr_fce)
      msg_error_f("number of independent AOs (%d) is not equal to"
       " number of spherical functions (%d)",1,d->n_ibfce,d->n_sphr_fce);
    n_blocks = (d->n_ibfce%n_vals ? d->n_ibfce/n_vals+1 : d->n_ibfce/n_vals);
    if (!d->ovrl)
      d->ovrl = mat_falloc(d->n_ibfce,d->n_ibfce);
    /* read data blocks */
    for (i=0; i<n_blocks; i++) {
      n_cols = ((i+1)<n_blocks ? n_vals :
       (d->n_ibfce%n_vals ? d->n_ibfce%n_vals : n_vals));
      i1 = 0;
      /* header */
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      line = str_free(line);
      /* atoms */
      for (j=0; j<d->n_atoms; j++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of cp2k log file",1);
        line = str_free(line);
        /* atom orbitals*/
        for (k=0; k<d->kind[d->atom[j].kind].n_sphr_fce; k++) {
          line = str_read_line_new(f);
          if (!line)
            msg_error("unexpected end of cp2k log file",1);
          /* overlap-matrix elements */
          i2 = n_vals*i;
          for (l=0; l<n_cols; l++)
            if (!sscanf(line+25+l*13,"%lf",&(d->ovrl[i1][i2++]))) 
              msg_error("invalid format of overlap matrix data block",1);
          line = str_free(line);
          i1++;
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
