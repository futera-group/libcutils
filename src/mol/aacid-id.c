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
#include "mol/aacid.h"

/* -------------------------------------------------------------------------- */

/* Convert one-letter amino-acid name to internal ID
   
   sym - the one-letter symbol */
short aacid_id_1(char sym) {
  switch (sym) {
    case 'A': return(AA_ALA);
    case 'R': return(AA_ARG);
    case 'N': return(AA_ASN);
    case 'D': return(AA_ASP);
    case 'C': return(AA_CYS);
    case 'E': return(AA_GLU);
    case 'Q': return(AA_GLN);
    case 'G': return(AA_GLY);
    case 'H': return(AA_HIS);
    case 'I': return(AA_ILE);
    case 'L': return(AA_LEU);
    case 'K': return(AA_LYS);
    case 'M': return(AA_MET);
    case 'F': return(AA_PHE);
    case 'P': return(AA_PRO);
    case 'S': return(AA_SER);
    case 'T': return(AA_THR);
    case 'W': return(AA_TRP);
    case 'Y': return(AA_TYR);
    case 'V': return(AA_VAL);
    }
  msg_error_f("invalid amino-acid symbol '%c'",1,sym);
  return(0);
  }

/* Convert three-letter amino-acid name to internal ID
   
   sym - the one-letter symbol */
short aacid_id_3(char *sym) {
  if      (str_compare_nc(sym,"ala")) return(AA_ALA);
  else if (str_compare_nc(sym,"arg")) return(AA_ARG);
  else if (str_compare_nc(sym,"asn")) return(AA_ASN);
  else if (str_compare_nc(sym,"asp")) return(AA_ASP);
  else if (str_compare_nc(sym,"cys")) return(AA_CYS);
  else if (str_compare_nc(sym,"glu")) return(AA_GLU);
  else if (str_compare_nc(sym,"gln")) return(AA_GLN);
  else if (str_compare_nc(sym,"gly")) return(AA_GLY);
  else if (str_compare_nc(sym,"his")) return(AA_HIS);
  else if (str_compare_nc(sym,"ile")) return(AA_ILE);
  else if (str_compare_nc(sym,"leu")) return(AA_LEU);
  else if (str_compare_nc(sym,"lys")) return(AA_LYS);
  else if (str_compare_nc(sym,"met")) return(AA_MET);
  else if (str_compare_nc(sym,"phe")) return(AA_PHE);
  else if (str_compare_nc(sym,"pro")) return(AA_PRO);
  else if (str_compare_nc(sym,"ser")) return(AA_SER);
  else if (str_compare_nc(sym,"thr")) return(AA_THR);
  else if (str_compare_nc(sym,"trp")) return(AA_TRP);
  else if (str_compare_nc(sym,"tyr")) return(AA_TYR);
  else if (str_compare_nc(sym,"val")) return(AA_VAL);
  else  
    msg_error_f("invalid amino-acid symbol '%s'",1,sym);
  return(0);
  }

/* Convert full amino-acid name to internal ID
   
   name - name of the amino acid */
short aacid_id_full(char *name) {
  if      (str_compare_nc(name,"alanine"      )) return(AA_ALA);
  else if (str_compare_nc(name,"arginine"     )) return(AA_ARG);
  else if (str_compare_nc(name,"asparagine"   )) return(AA_ASN);
  else if (str_compare_nc(name,"aspartic acid")) return(AA_ASP);
  else if (str_compare_nc(name,"cysteine"     )) return(AA_CYS);
  else if (str_compare_nc(name,"glutamic acid")) return(AA_GLU);
  else if (str_compare_nc(name,"glutamine"    )) return(AA_GLN);
  else if (str_compare_nc(name,"glycine"      )) return(AA_GLY);
  else if (str_compare_nc(name,"histidine"    )) return(AA_HIS);
  else if (str_compare_nc(name,"isoleucine"   )) return(AA_ILE);
  else if (str_compare_nc(name,"leucine"      )) return(AA_LEU);
  else if (str_compare_nc(name,"lysine"       )) return(AA_LYS);
  else if (str_compare_nc(name,"methionine"   )) return(AA_MET);
  else if (str_compare_nc(name,"phenylalanine")) return(AA_PHE);
  else if (str_compare_nc(name,"proline"      )) return(AA_PRO);
  else if (str_compare_nc(name,"serine"       )) return(AA_SER);
  else if (str_compare_nc(name,"threonine"    )) return(AA_THR);
  else if (str_compare_nc(name,"tryptophan"   )) return(AA_TRP);
  else if (str_compare_nc(name,"tyrosine"     )) return(AA_TYR);
  else if (str_compare_nc(name,"valine"       )) return(AA_VAL);
  else  
    msg_error_f("invalid amino-acid name '%s'",1,name);
  return(0);
  }

/* -------------------------------------------------------------------------- */
