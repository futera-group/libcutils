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
#include <cmn/vector.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* assing IOp value to internal code array

   c - pointer to IOp code array
   o - overlay number 
   p - internal option number
   v - internal option value */
void gauss_iop_assign(long *c, unsigned o, unsigned p, int v) {
  switch (o) {
    case  1:
      switch (p) {
        case  38: c[GAUSS_IOP_01_038]=v; break;
        case 172: c[GAUSS_IOP_01_172]=v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  2:
      switch (p) {
        case  12: c[GAUSS_IOP_02_012] = v; break;
        case  17: c[GAUSS_IOP_02_017] = v; break;
        case  18: c[GAUSS_IOP_02_018] = v; break;
        case  40: c[GAUSS_IOP_02_040] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  3:
      switch (p) {
        case   5: c[GAUSS_IOP_03_005] = v; break;
        case   6: c[GAUSS_IOP_03_006] = v; break;
        case   7: c[GAUSS_IOP_03_007] = v; break;
        case  11: c[GAUSS_IOP_03_011] = v; break;
        case  16: c[GAUSS_IOP_03_016] = v; break;
        case  17: c[GAUSS_IOP_03_017] = v; break;
        case  24: c[GAUSS_IOP_03_024] = v; break;
        case  25: c[GAUSS_IOP_03_025] = v; break;
        case  30: c[GAUSS_IOP_03_030] = v; break;
        case  41: c[GAUSS_IOP_03_041] = v; break;
        case  74: c[GAUSS_IOP_03_074] = v; break;
        case  82: c[GAUSS_IOP_03_082] = v; break;
        case 116: c[GAUSS_IOP_03_116] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  4:
      switch (p) {
        case   5: c[GAUSS_IOP_04_005] = v; break;
        case  35: c[GAUSS_IOP_04_035] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  5:
      switch (p) {
        case   5: c[GAUSS_IOP_05_005] = v; break;
        case   7: c[GAUSS_IOP_05_007] = v; break;
        case   8: c[GAUSS_IOP_05_008] = v; break;
        case  13: c[GAUSS_IOP_05_013] = v; break;
        case  32: c[GAUSS_IOP_05_032] = v; break;
        case  35: c[GAUSS_IOP_05_035] = v; break;
        case  38: c[GAUSS_IOP_05_038] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  6:
      switch (p) {
        case   7: c[GAUSS_IOP_06_007] = v; break;
        case   8: c[GAUSS_IOP_06_008] = v; break;
        case   9: c[GAUSS_IOP_06_009] = v; break;
        case  10: c[GAUSS_IOP_06_010] = v; break;
        case  28: c[GAUSS_IOP_06_028] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  8:
      switch (p) {
        case   6: c[GAUSS_IOP_08_006] = v; break;
        case   9: c[GAUSS_IOP_08_009] = v; break;
        case  10: c[GAUSS_IOP_08_010] = v; break;
        case  68: c[GAUSS_IOP_08_068] = v; break;
        case 108: c[GAUSS_IOP_08_108] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case  9:
      switch (p) {
        case   5: c[GAUSS_IOP_09_005] = v; break;
        case  23: c[GAUSS_IOP_09_023] = v; break;
        case  41: c[GAUSS_IOP_09_041] = v; break;
        case  42: c[GAUSS_IOP_09_042] = v; break;
        case  44: c[GAUSS_IOP_09_044] = v; break;
        case  49: c[GAUSS_IOP_09_049] = v; break;
        case  68: c[GAUSS_IOP_09_068] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    case 99:
      switch (p) {
        case   5: c[GAUSS_IOP_99_005] = v; break;
        case   9: c[GAUSS_IOP_99_009] = v; break;
        default:
          msg_warn_f("undefined IOp %d in overlay %d",p,o);
        }
      break;
    default:
      msg_warn_f("undefined overlay %d for IOp %d",o,p);
    }
  }

/* set IOp codes from one line specification
 
   p - pointer to IOp code array
   s - line from gaussian output file */
void gauss_iop_set(long *p, char *s) {
  unsigned i,ovr,iop,iop_n,pair_n;
  char **iop_s,**pair_s;
  long val;
  if (sscanf(s,"%u",&ovr)!=1)
    msg_error("cannot read overlay number from IPs line",1);
  str_trim_mark(s,'/');
  iop_s = str_split(s,',',&iop_n);
  for (i=0; i<iop_n; i++) {
    pair_s = str_split(iop_s[i],'=',&pair_n);
    if (pair_n!=2 ||
        sscanf(pair_s[0],"%u",&iop)!=1 ||
        sscanf(pair_s[1],"%ld",&val)!=1)
      msg_error("invalid format of IOp specification",1);
    gauss_iop_assign(p,ovr,iop,val);
    vec_sfree(pair_s,pair_n);
    }
  vec_sfree(iop_s,iop_n);
  }

/* -------------------------------------------------------------------------- */
