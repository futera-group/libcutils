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

#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert tree mark to internal code
 
   t - the tree mark from file */
short mol_apc_tree_id(char *t) {
  if (str_compare(t,"M"))
    return(TREE_MARK_M);
  else if (str_compare(t,"E"))
    return(TREE_MARK_E);
  else if (str_compare(t,"S"))
    return(TREE_MARK_S);
  else if (str_compare(t,"B"))
    return(TREE_MARK_B);
  else if (str_compare(t,"3"))
    return(TREE_MARK_3);
  else if (str_compare(t,"4"))
    return(TREE_MARK_4);
  else if (str_compare(t,"5"))
    return(TREE_MARK_5);
  else if (str_compare(t,"6"))
    return(TREE_MARK_6);
  else
    msg_error_f("unknown tree mark \"%s\"",1,t);
  return(0);
  }

/* convert tree internal code to tree mark
 
   t - the tree internal code
   m - the tree mark */
void mol_apc_tree_mark(short t, char *m) {
  switch (t) {
    case TREE_MARK_M: sprintf(m,"M"); break;
    case TREE_MARK_E: sprintf(m,"E"); break;
    case TREE_MARK_S: sprintf(m,"S"); break;
    case TREE_MARK_B: sprintf(m,"B"); break;
    case TREE_MARK_3: sprintf(m,"3"); break;
    case TREE_MARK_4: sprintf(m,"4"); break;
    case TREE_MARK_5: sprintf(m,"5"); break;
    case TREE_MARK_6: sprintf(m,"6"); break;
    default:
      msg_error("invalid internal tree mark identificator",1);
      break;
    }
  }

/* -------------------------------------------------------------------------- */
