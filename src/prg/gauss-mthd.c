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

#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* analyze control line with computational options and return method ID
 
   s - line from gaussian output file */
short gauss_mthd_id(char *s) {
  char **cmd_s,**opt_s,key[80];
  short mthd_id = GAUSS_MTHD_UNK;
  unsigned i,j,cmd_n,opt_n;
  cmd_s = str_split(s,' ',&cmd_n);
  for (i=0; i<cmd_n; i++) {
    if (cmd_s[i][0]!='#') {
      str_lowcase(cmd_s[i]);
      for (j=0; j<str_length(cmd_s[i]) && 
           cmd_s[i][j]!='=' && cmd_s[i][j]!='(' && cmd_s[i][j]!='/'; j++)
        key[j] = cmd_s[i][j];
      key[j] = '\0';
      /* post-HF methods */
      if (mthd_id==GAUSS_MTHD_HF) {
        if (str_compare(key,"cis")) {
          mthd_id = GAUSS_MTHD_CIS;
          str_ltrim_mark(cmd_s[i],'(');
          str_rtrim_mark(cmd_s[i],')');
          opt_s = str_split(cmd_s[i],',',&opt_n);
          for (j=0; j<opt_n; j++)
            if (str_compare(opt_s[j],"d")) {
              mthd_id = GAUSS_MTHD_CISd;
              }
          vec_sfree(opt_s,opt_n);
          }
        }
      /* other methods */
      else {
        if (str_compare(key,"hf"))
          mthd_id = GAUSS_MTHD_HF;
        else if (str_compare(key,"eomccsd"))
          mthd_id = GAUSS_MTHD_EOMCCSD;
        else if (str_compare(key,"b3lyp"))
          mthd_id = GAUSS_MTHD_DFT;
        }
      }
    }
  vec_sfree(cmd_s,cmd_n);
  return(mthd_id);
  }

/* -------------------------------------------------------------------------- */
