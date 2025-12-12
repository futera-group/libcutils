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
#include <stdlib.h>
#include <cmn/string.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* return string with bibtex entry name generated from the saved data
 
   e - the to bibtex entry */
char* bibtex_entry_name_new(struct bib_entry *e) {
  char *label = NULL,*s;
  if (e) {
    /* first author + year */
    if (e->author && e->n_authors) {
      if (e->year) {
        label = str_new(str_length(e->author[0][0])+str_unum_len(e->year));
        sprintf(label,"%s%u",e->author[0][0],e->year);
        }
      else
        label = str_copy_new(e->author[0][0]);
      }
    /* first editor + year */
    else if (e->editor && e->n_editors) {
      if (e->year) {
        label = str_new(str_length(e->editor[0][0])+str_unum_len(e->year));
        sprintf(label,"%s%u",e->editor[0][0],e->year);
        }
      else
        label = str_copy_new(e->editor[0][0]);
      }
    /* record type + year */
    else {
      s = bibtex_type_entry_name(e->type);
      str_Lowcase(s);
      if (e->year) {
        label = str_new(str_length(s)+str_unum_len(e->year));
        sprintf(label,"%s%u",s,e->year);
        }
      else
        label = str_copy_new(s);
      }
    }
  return(label);
  }

/* -------------------------------------------------------------------------- */
