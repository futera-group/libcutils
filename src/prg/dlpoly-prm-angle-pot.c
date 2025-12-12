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
#include <cmn/string.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of angle potential
 
   s - the potential symbol */
short dlpoly_prm_angle_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"harm")) return(DLPOLY_ANGL_HM);
  if (str_compare(s,"-hrm")) return(DLPOLY_ANGL_HM);
  if (str_compare(s,"quar")) return(DLPOLY_ANGL_QR);
  if (str_compare(s,"-qur")) return(DLPOLY_ANGL_QR);
  if (str_compare(s,"thrm")) return(DLPOLY_ANGL_TH);
  if (str_compare(s,"-thm")) return(DLPOLY_ANGL_TH);
  if (str_compare(s,"shrm")) return(DLPOLY_ANGL_SH);
  if (str_compare(s,"-shm")) return(DLPOLY_ANGL_SH);
  if (str_compare(s,"bvs1")) return(DLPOLY_ANGL_SV);
  if (str_compare(s,"-bv1")) return(DLPOLY_ANGL_SV);
  if (str_compare(s,"bvs2")) return(DLPOLY_ANGL_TV);
  if (str_compare(s,"-bv2")) return(DLPOLY_ANGL_TV);
  if (str_compare(s,"hcos")) return(DLPOLY_ANGL_HC);
  if (str_compare(s,"-hcs")) return(DLPOLY_ANGL_HC);
  if (str_compare(s,"cos"))  return(DLPOLY_ANGL_CS);
  if (str_compare(s,"-cos")) return(DLPOLY_ANGL_CS);
  if (str_compare(s,"mmsb")) return(DLPOLY_ANGL_MM);
  if (str_compare(s,"-msb")) return(DLPOLY_ANGL_MM);
  if (str_compare(s,"stst")) return(DLPOLY_ANGL_SS);
  if (str_compare(s,"-sts")) return(DLPOLY_ANGL_SS);
  if (str_compare(s,"stbe")) return(DLPOLY_ANGL_SB);
  if (str_compare(s,"-stb")) return(DLPOLY_ANGL_SB);
  if (str_compare(s,"cmps")) return(DLPOLY_ANGL_CM);
  if (str_compare(s,"-cmp")) return(DLPOLY_ANGL_CM);
  if (str_compare(s,"amoe")) return(DLPOLY_ANGL_AM);
  if (str_compare(s,"-amo")) return(DLPOLY_ANGL_AM);
  if (str_compare(s,"kky"))  return(DLPOLY_ANGL_KK);
  if (str_compare(s,"-kky")) return(DLPOLY_ANGL_KK);
  return(id);
  }

/* convert internal ID of angle potential to its symbol
 
   id - the internal ID
   m  - minus form of the keyword
   s  - storage for the symbol */
char *dlpoly_prm_angle_pot_sym(short id, short m) {
  static char s[80];
  if (m) {
    /* full form keywords */
    switch (id) {
      case DLPOLY_ANGL_HM: sprintf(s,"%s","-hrm"); break;
      case DLPOLY_ANGL_QR: sprintf(s,"%s","-qur"); break;
      case DLPOLY_ANGL_TH: sprintf(s,"%s","-thm"); break;
      case DLPOLY_ANGL_SH: sprintf(s,"%s","-shm"); break;
      case DLPOLY_ANGL_SV: sprintf(s,"%s","-bv1"); break;
      case DLPOLY_ANGL_TV: sprintf(s,"%s","-bv2"); break;
      case DLPOLY_ANGL_HC: sprintf(s,"%s","-hcs"); break;
      case DLPOLY_ANGL_CS: sprintf(s,"%s","-cos"); break;
      case DLPOLY_ANGL_MM: sprintf(s,"%s","-msb"); break;
      case DLPOLY_ANGL_SS: sprintf(s,"%s","-sts"); break;
      case DLPOLY_ANGL_SB: sprintf(s,"%s","-stb"); break;
      case DLPOLY_ANGL_CM: sprintf(s,"%s","-cmp"); break;
      case DLPOLY_ANGL_AM: sprintf(s,"%s","-amo"); break;
      case DLPOLY_ANGL_KK: sprintf(s,"%s","-kky"); break;
      }
    }
  else {
    /* minus form keywords */
    switch (id) {
      case DLPOLY_ANGL_HM: sprintf(s,"%s","harm"); break;
      case DLPOLY_ANGL_QR: sprintf(s,"%s","quar"); break;
      case DLPOLY_ANGL_TH: sprintf(s,"%s","thrm"); break;
      case DLPOLY_ANGL_SH: sprintf(s,"%s","shrm"); break;
      case DLPOLY_ANGL_SV: sprintf(s,"%s","bvs1"); break;
      case DLPOLY_ANGL_TV: sprintf(s,"%s","bvs2"); break;
      case DLPOLY_ANGL_HC: sprintf(s,"%s","hcos"); break;
      case DLPOLY_ANGL_CS: sprintf(s,"%s","cos");  break;
      case DLPOLY_ANGL_MM: sprintf(s,"%s","mmsb"); break;
      case DLPOLY_ANGL_SS: sprintf(s,"%s","stst"); break;
      case DLPOLY_ANGL_SB: sprintf(s,"%s","stbe"); break;
      case DLPOLY_ANGL_CM: sprintf(s,"%s","cmps"); break;
      case DLPOLY_ANGL_AM: sprintf(s,"%s","amoe"); break;
      case DLPOLY_ANGL_KK: sprintf(s,"%s","kky");  break;
      }
    }
  return(s);
  }

/* return number of variables of specific angle potential
 
   id - the potential specifier */
unsigned dlpoly_prm_angle_pot_nvar(short id) {
  unsigned n = 0;
  switch (id) {
    case DLPOLY_ANGL_HM: n = 2; break;
    case DLPOLY_ANGL_QR: n = 4; break;
    case DLPOLY_ANGL_TH: n = 3; break;
    case DLPOLY_ANGL_SH: n = 4; break;
    case DLPOLY_ANGL_SV: n = 4; break;
    case DLPOLY_ANGL_TV: n = 4; break;
    case DLPOLY_ANGL_HC: n = 2; break;
    case DLPOLY_ANGL_CS: n = 3; break;
    case DLPOLY_ANGL_MM: n = 4; break;
    case DLPOLY_ANGL_SS: n = 3; break;
    case DLPOLY_ANGL_SB: n = 3; break;
    case DLPOLY_ANGL_CM: n = 4; break;
    case DLPOLY_ANGL_AM: n = 2; break;
    case DLPOLY_ANGL_KK: n = 4; break;
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
