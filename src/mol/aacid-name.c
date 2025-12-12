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
#include "mol/aacid.h"

/* -------------------------------------------------------------------------- */

/* Return full name of specified amino acid
   
   id - internal amino-acid ID */
char* aacid_name_full(short id) {
  static char name[80];
  name[0] = '\0';
  switch (id) {
    case AA_ALA: sprintf(name,"%s","alanine");       break;
    case AA_ARG: sprintf(name,"%s","arginine");      break;
    case AA_ASN: sprintf(name,"%s","asparagine");    break;
    case AA_ASP: sprintf(name,"%s","aspartic acid"); break;
    case AA_CYS: sprintf(name,"%s","cysteine");      break;
    case AA_GLU: sprintf(name,"%s","glutamic acid"); break;
    case AA_GLN: sprintf(name,"%s","glutamine");     break;
    case AA_GLY: sprintf(name,"%s","glycine");       break;
    case AA_HIS: sprintf(name,"%s","histidine");     break;
    case AA_ILE: sprintf(name,"%s","isoleucine");    break;
    case AA_LEU: sprintf(name,"%s","leucine");       break;
    case AA_LYS: sprintf(name,"%s","lysine");        break;
    case AA_MET: sprintf(name,"%s","methionine");    break;
    case AA_PHE: sprintf(name,"%s","phenylalanine"); break;
    case AA_PRO: sprintf(name,"%s","proline");       break;
    case AA_SER: sprintf(name,"%s","serine");        break;
    case AA_THR: sprintf(name,"%s","threonine");     break;
    case AA_TRP: sprintf(name,"%s","tryptophan");    break;
    case AA_TYR: sprintf(name,"%s","tyrosine");      break;
    case AA_VAL: sprintf(name,"%s","valine");        break;
    }
  return(name);
  }

/* Return 3-letter name of specified amino acid
   
   id - internal amino-acid ID */
char* aacid_name_3(short id) {
  static char name[5];
  name[0] = '\0';
  switch (id) {
    case AA_ALA: sprintf(name,"%s","Ala"); break;
    case AA_ARG: sprintf(name,"%s","Arg"); break;
    case AA_ASN: sprintf(name,"%s","Asn"); break;
    case AA_ASP: sprintf(name,"%s","Asp"); break;
    case AA_CYS: sprintf(name,"%s","Cys"); break;
    case AA_GLU: sprintf(name,"%s","Glu"); break;
    case AA_GLN: sprintf(name,"%s","Gln"); break;
    case AA_GLY: sprintf(name,"%s","Gly"); break;
    case AA_HIS: sprintf(name,"%s","His"); break;
    case AA_ILE: sprintf(name,"%s","Ile"); break;
    case AA_LEU: sprintf(name,"%s","Leu"); break;
    case AA_LYS: sprintf(name,"%s","Lys"); break;
    case AA_MET: sprintf(name,"%s","Met"); break;
    case AA_PHE: sprintf(name,"%s","Phe"); break;
    case AA_PRO: sprintf(name,"%s","Pro"); break;
    case AA_SER: sprintf(name,"%s","Ser"); break;
    case AA_THR: sprintf(name,"%s","Thr"); break;
    case AA_TRP: sprintf(name,"%s","Trp"); break;
    case AA_TYR: sprintf(name,"%s","Tyr"); break;
    case AA_VAL: sprintf(name,"%s","Val"); break;
    }
  return(name);
  }

/* Return 1-letter symbol of specified amino acid
   
   id - internal amino-acid ID */
char* aacid_name_1(short id) {
  static char name[5];
  name[0] = '\0';
  switch (id) {
    case AA_ALA: sprintf(name,"%s","A"); break;
    case AA_ARG: sprintf(name,"%s","R"); break;
    case AA_ASN: sprintf(name,"%s","N"); break;
    case AA_ASP: sprintf(name,"%s","D"); break;
    case AA_CYS: sprintf(name,"%s","C"); break;
    case AA_GLU: sprintf(name,"%s","E"); break;
    case AA_GLN: sprintf(name,"%s","Q"); break;
    case AA_GLY: sprintf(name,"%s","G"); break;
    case AA_HIS: sprintf(name,"%s","H"); break;
    case AA_ILE: sprintf(name,"%s","I"); break;
    case AA_LEU: sprintf(name,"%s","L"); break;
    case AA_LYS: sprintf(name,"%s","K"); break;
    case AA_MET: sprintf(name,"%s","M"); break;
    case AA_PHE: sprintf(name,"%s","F"); break;
    case AA_PRO: sprintf(name,"%s","P"); break;
    case AA_SER: sprintf(name,"%s","S"); break;
    case AA_THR: sprintf(name,"%s","T"); break;
    case AA_TRP: sprintf(name,"%s","W"); break;
    case AA_TYR: sprintf(name,"%s","Y"); break;
    case AA_VAL: sprintf(name,"%s","V"); break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
