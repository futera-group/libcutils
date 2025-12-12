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
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of pair potential
 
   s - the potential symbol */
short lammps_pot_pair_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"none"))                           return(LAMMPS_PAIR_NONE);
  if (str_compare(s,"zero"))                           return(LAMMPS_PAIR_ZERO);
  if (str_compare(s,"hybrid"))                         return(LAMMPS_PAIR_HYBR);
  if (str_compare(s,"hybrid/overlay"))                 return(LAMMPS_PAIR_HBOV);
  if (str_compare(s,"adp"))                            return(LAMMPS_PAIR_ADPT);
  if (str_compare(s,"airebo"))                         return(LAMMPS_PAIR_AIRB);
  if (str_compare(s,"airebo/morse"))                   return(LAMMPS_PAIR_AIRM);
  if (str_compare(s,"beck"))                           return(LAMMPS_PAIR_BECK);
  if (str_compare(s,"body"))                           return(LAMMPS_PAIR_BODY);
  if (str_compare(s,"bop"))                            return(LAMMPS_PAIR_BOPT);
  if (str_compare(s,"born"))                           return(LAMMPS_PAIR_BORN);
  if (str_compare(s,"born/coul/long"))                 return(LAMMPS_PAIR_BRCL);
  if (str_compare(s,"born/coul/long/cs"))              return(LAMMPS_PAIR_BRCC);
  if (str_compare(s,"born/coul/msm"))                  return(LAMMPS_PAIR_BRCM);
  if (str_compare(s,"born/coul/wolf"))                 return(LAMMPS_PAIR_BRCW);
  if (str_compare(s,"brownian"))                       return(LAMMPS_PAIR_BRWN);
  if (str_compare(s,"brownian/poly"))                  return(LAMMPS_PAIR_BRPL);
  if (str_compare(s,"buck"))                           return(LAMMPS_PAIR_BUCK);
  if (str_compare(s,"buck/coul/cut"))                  return(LAMMPS_PAIR_BCCT);
  if (str_compare(s,"buck/coul/long"))                 return(LAMMPS_PAIR_BCCL);
  if (str_compare(s,"buck/coul/long/cs"))              return(LAMMPS_PAIR_BCCC);
  if (str_compare(s,"buck/coul/msm"))                  return(LAMMPS_PAIR_BCCM);
  if (str_compare(s,"buck/long/coul/long"))            return(LAMMPS_PAIR_BCLC);
  if (str_compare(s,"colloid"))                        return(LAMMPS_PAIR_CLLD);
  if (str_compare(s,"comb"))                           return(LAMMPS_PAIR_COMB);
  if (str_compare(s,"comb3"))                          return(LAMMPS_PAIR_CMB3);
  if (str_compare(s,"coul/cut"))                       return(LAMMPS_PAIR_CLCT);
  if (str_compare(s,"coul/debye"))                     return(LAMMPS_PAIR_CLDB);
  if (str_compare(s,"coul/dsf"))                       return(LAMMPS_PAIR_CLDS);
  if (str_compare(s,"coul/long"))                      return(LAMMPS_PAIR_CLLG);
  if (str_compare(s,"coul/long/cs"))                   return(LAMMPS_PAIR_CLLC);
  if (str_compare(s,"coul/msm"))                       return(LAMMPS_PAIR_CLMS);
  if (str_compare(s,"coul/streitz"))                   return(LAMMPS_PAIR_CLSM);
  if (str_compare(s,"could/wolf"))                     return(LAMMPS_PAIR_CLWF);
  if (str_compare(s,"dpd"))                            return(LAMMPS_PAIR_DPDP);
  if (str_compare(s,"dpd/tstat"))                      return(LAMMPS_PAIR_DPDT);
  if (str_compare(s,"dsmc"))                           return(LAMMPS_PAIR_DSMC);
  if (str_compare(s,"eam"))                            return(LAMMPS_PAIR_EAMP);
  if (str_compare(s,"eam/alloy"))                      return(LAMMPS_PAIR_EAMA);
  if (str_compare(s,"eam/fs"))                         return(LAMMPS_PAIR_EAMF);
  if (str_compare(s,"eim"))                            return(LAMMPS_PAIR_EIMP);
  if (str_compare(s,"gauss"))                          return(LAMMPS_PAIR_GAUS);
  if (str_compare(s,"gayberne"))                       return(LAMMPS_PAIR_GBEP);
  if (str_compare(s,"gran/hertz/history"))             return(LAMMPS_PAIR_GRHR);
  if (str_compare(s,"gran/hooke"))                     return(LAMMPS_PAIR_GRHK);
  if (str_compare(s,"gran/hooke/history"))             return(LAMMPS_PAIR_GRHH);
  if (str_compare(s,"hbond/dreiding/lj"))              return(LAMMPS_PAIR_HBDL);
  if (str_compare(s,"hbond/dreiding/morse"))           return(LAMMPS_PAIR_HBDM);
  if (str_compare(s,"kim"))                            return(LAMMPS_PAIR_KIMP);
  if (str_compare(s,"lcbop"))                          return(LAMMPS_PAIR_CBOP);
  if (str_compare(s,"line/lj"))                        return(LAMMPS_PAIR_LJPT);
  if (str_compare(s,"lj/charmm/coul/charmm"))          return(LAMMPS_PAIR_CHCC);
  if (str_compare(s,"lj/charmm/coul/charmm/implicit")) return(LAMMPS_PAIR_CHCI);
  if (str_compare(s,"lj/charmm/coul/long"))            return(LAMMPS_PAIR_CHCL);
  if (str_compare(s,"lj/charmm/coul/msm"))             return(LAMMPS_PAIR_CHCM);
  if (str_compare(s,"lj/class2"))                      return(LAMMPS_PAIR_CLS2);
  if (str_compare(s,"lj/class2/coul/cut"))             return(LAMMPS_PAIR_C2CC);
  if (str_compare(s,"lj/class2/coul/long"))            return(LAMMPS_PAIR_C2LC);
  if (str_compare(s,"lj/cubic"))                       return(LAMMPS_PAIR_LJCB);
  if (str_compare(s,"lj/cut"))                         return(LAMMPS_PAIR_LJCT);
  if (str_compare(s,"lj/cut/coul/cut"))                return(LAMMPS_PAIR_LJCC);
  if (str_compare(s,"lj/cut/coul/debye"))              return(LAMMPS_PAIR_LJCD);
  if (str_compare(s,"lj/cut/coul/dsf"))                return(LAMMPS_PAIR_LJCS);
  if (str_compare(s,"lj/cut/coul/long"))               return(LAMMPS_PAIR_LJCL);
  if (str_compare(s,"lj/cut/coul/long/cs"))            return(LAMMPS_PAIR_LJCR);
  if (str_compare(s,"lj/cut/coul/msm"))                return(LAMMPS_PAIR_LJCM);
  if (str_compare(s,"lj/cut/dipole/cut"))              return(LAMMPS_PAIR_LJDC);
  if (str_compare(s,"lj/cut/dipole/long"))             return(LAMMPS_PAIR_LJDL);
  if (str_compare(s,"lj/cut/tip4p/cut"))               return(LAMMPS_PAIR_LJT4);
  if (str_compare(s,"lj/cut/tip4p/long"))              return(LAMMPS_PAIR_LJL4);
  if (str_compare(s,"lj/expand"))                      return(LAMMPS_PAIR_LJVS);
  if (str_compare(s,"lj/gromacs"))                     return(LAMMPS_PAIR_GXLJ);
  if (str_compare(s,"lj/gromacs/coul/gromacs"))        return(LAMMPS_PAIR_GXCL);
  if (str_compare(s,"lj/long/coul/long"))              return(LAMMPS_PAIR_LLCL);
  if (str_compare(s,"lj/long/dipole/long"))            return(LAMMPS_PAIR_LLDL);
  if (str_compare(s,"lj/long/tip4p/long"))             return(LAMMPS_PAIR_LLT4);
  if (str_compare(s,"lj/smooth"))                      return(LAMMPS_PAIR_SMLJ);
  if (str_compare(s,"lj/smooth/linear"))               return(LAMMPS_PAIR_SMLN);
  if (str_compare(s,"lj96/cut"))                       return(LAMMPS_PAIR_LJ96);
  if (str_compare(s,"lubricate"))                      return(LAMMPS_PAIR_LBPT);
  if (str_compare(s,"lubricate/poly"))                 return(LAMMPS_PAIR_LBPL);
  if (str_compare(s,"lubricateU"))                     return(LAMMPS_PAIR_LBFD);
  if (str_compare(s,"lubricateU/poly"))                return(LAMMPS_PAIR_LBFP);
  if (str_compare(s,"meam"))                           return(LAMMPS_PAIR_MEAM);
  if (str_compare(s,"mie/cut"))                        return(LAMMPS_PAIR_MIEP);
  if (str_compare(s,"morse"))                          return(LAMMPS_PAIR_MORS);
  if (str_compare(s,"nb3b/harmonic"))                  return(LAMMPS_PAIR_NB3B);
  if (str_compare(s,"nm/cut"))                         return(LAMMPS_PAIR_NMPT);
  if (str_compare(s,"nm/cut/coul/cut"))                return(LAMMPS_PAIR_NMCC);
  if (str_compare(s,"nm/cut/coul/long"))               return(LAMMPS_PAIR_NMCL);
  if (str_compare(s,"peri/eps"))                       return(LAMMPS_PAIR_PREP);
  if (str_compare(s,"peri/lps"))                       return(LAMMPS_PAIR_PRLP);
  if (str_compare(s,"peri/pmb"))                       return(LAMMPS_PAIR_PRPM);
  if (str_compare(s,"peri/ves"))                       return(LAMMPS_PAIR_PRVS);
  if (str_compare(s,"polymorphic"))                    return(LAMMPS_PAIR_PLMR);
  if (str_compare(s,"reax"))                           return(LAMMPS_PAIR_REAX);
  if (str_compare(s,"rebo"))                           return(LAMMPS_PAIR_REBO);
  if (str_compare(s,"resquared"))                      return(LAMMPS_PAIR_RESQ);
  if (str_compare(s,"snap"))                           return(LAMMPS_PAIR_SNAP);
  if (str_compare(s,"soft"))                           return(LAMMPS_PAIR_SOFT);
  if (str_compare(s,"sw"))                             return(LAMMPS_PAIR_SWPT);
  if (str_compare(s,"table"))                          return(LAMMPS_PAIR_TABL);
  if (str_compare(s,"tersoff"))                        return(LAMMPS_PAIR_TERS);
  if (str_compare(s,"tersoff/mod"))                    return(LAMMPS_PAIR_TERM);
  if (str_compare(s,"tersoff/zbi"))                    return(LAMMPS_PAIR_TERZ);
  if (str_compare(s,"tip4p/cut"))                      return(LAMMPS_PAIR_TP4C);
  if (str_compare(s,"tip4p/long"))                     return(LAMMPS_PAIR_TP4L);
  if (str_compare(s,"tri/lj"))                         return(LAMMPS_PAIR_TRIL);
  if (str_compare(s,"vashishta"))                      return(LAMMPS_PAIR_VASH);
  if (str_compare(s,"yukawa"))                         return(LAMMPS_PAIR_YUKW);
  if (str_compare(s,"yukawa/colloid"))                 return(LAMMPS_PAIR_YUKC);
  if (str_compare(s,"zbl"))                            return(LAMMPS_PAIR_ZBLP);
  msg_error_f("unknown pair potential type \"%s\"",1,s);
  return(id);
  }

/* return symbol of pair potential
 
   t - internal ID of the potential
   s - the potential symbol */
char* lammps_pot_pair_sym(short t) {
  static char s[80];
  s[0]='\0';
  switch (t) {
    case LAMMPS_PAIR_NONE: sprintf(s,"none");                           break;
    case LAMMPS_PAIR_ZERO: sprintf(s,"zero");                           break;
    case LAMMPS_PAIR_HYBR: sprintf(s,"hybrid");                         break;
    case LAMMPS_PAIR_HBOV: sprintf(s,"hybrid/overlay");                 break;
    case LAMMPS_PAIR_ADPT: sprintf(s,"adp");                            break;
    case LAMMPS_PAIR_AIRB: sprintf(s,"airebo");                         break;
    case LAMMPS_PAIR_AIRM: sprintf(s,"airebo/morse");                   break;
    case LAMMPS_PAIR_BECK: sprintf(s,"beck");                           break;
    case LAMMPS_PAIR_BODY: sprintf(s,"body");                           break;
    case LAMMPS_PAIR_BOPT: sprintf(s,"bop");                            break;
    case LAMMPS_PAIR_BORN: sprintf(s,"born");                           break;
    case LAMMPS_PAIR_BRCL: sprintf(s,"born/coul/long");                 break;
    case LAMMPS_PAIR_BRCC: sprintf(s,"born/coul/long/cs");              break;
    case LAMMPS_PAIR_BRCM: sprintf(s,"born/coul/msm");                  break;
    case LAMMPS_PAIR_BRCW: sprintf(s,"born/coul/wolf");                 break;
    case LAMMPS_PAIR_BRWN: sprintf(s,"brownian");                       break;
    case LAMMPS_PAIR_BRPL: sprintf(s,"brownian/poly");                  break;
    case LAMMPS_PAIR_BUCK: sprintf(s,"buck");                           break;
    case LAMMPS_PAIR_BCCT: sprintf(s,"buck/coul/cut");                  break;
    case LAMMPS_PAIR_BCCL: sprintf(s,"buck/coul/long");                 break;
    case LAMMPS_PAIR_BCCC: sprintf(s,"buck/coul/long/cs");              break;
    case LAMMPS_PAIR_BCCM: sprintf(s,"buck/coul/msm");                  break;
    case LAMMPS_PAIR_BCLC: sprintf(s,"buck/long/coul/long");            break;
    case LAMMPS_PAIR_CLLD: sprintf(s,"colloid");                        break;
    case LAMMPS_PAIR_COMB: sprintf(s,"comb");                           break;
    case LAMMPS_PAIR_CMB3: sprintf(s,"comb3");                          break;
    case LAMMPS_PAIR_CLCT: sprintf(s,"coul/cut");                       break;
    case LAMMPS_PAIR_CLDB: sprintf(s,"coul/debye");                     break;
    case LAMMPS_PAIR_CLDS: sprintf(s,"coul/dsf");                       break;
    case LAMMPS_PAIR_CLLG: sprintf(s,"coul/long");                      break;
    case LAMMPS_PAIR_CLLC: sprintf(s,"coul/long/cs");                   break;
    case LAMMPS_PAIR_CLMS: sprintf(s,"coul/msm");                       break;
    case LAMMPS_PAIR_CLSM: sprintf(s,"coul/streitz");                   break;
    case LAMMPS_PAIR_CLWF: sprintf(s,"could/wolf");                     break;
    case LAMMPS_PAIR_DPDP: sprintf(s,"dpd");                            break;
    case LAMMPS_PAIR_DPDT: sprintf(s,"dpd/tstat");                      break;
    case LAMMPS_PAIR_DSMC: sprintf(s,"dsmc");                           break;
    case LAMMPS_PAIR_EAMP: sprintf(s,"eam");                            break;
    case LAMMPS_PAIR_EAMA: sprintf(s,"eam/alloy");                      break;
    case LAMMPS_PAIR_EAMF: sprintf(s,"eam/fs");                         break;
    case LAMMPS_PAIR_EIMP: sprintf(s,"eim");                            break;
    case LAMMPS_PAIR_GAUS: sprintf(s,"gauss");                          break;
    case LAMMPS_PAIR_GBEP: sprintf(s,"gayberne");                       break;
    case LAMMPS_PAIR_GRHR: sprintf(s,"gran/hertz/history");             break;
    case LAMMPS_PAIR_GRHK: sprintf(s,"gran/hooke");                     break;
    case LAMMPS_PAIR_GRHH: sprintf(s,"gran/hooke/history");             break;
    case LAMMPS_PAIR_HBDL: sprintf(s,"hbond/dreiding/lj");              break;
    case LAMMPS_PAIR_HBDM: sprintf(s,"hbond/dreiding/morse");           break;
    case LAMMPS_PAIR_KIMP: sprintf(s,"kim");                            break;
    case LAMMPS_PAIR_CBOP: sprintf(s,"lcbop");                          break;
    case LAMMPS_PAIR_LJPT: sprintf(s,"line/lj");                        break;
    case LAMMPS_PAIR_CHCC: sprintf(s,"lj/charmm/coul/charmm");          break;
    case LAMMPS_PAIR_CHCI: sprintf(s,"lj/charmm/coul/charmm/implicit"); break;
    case LAMMPS_PAIR_CHCL: sprintf(s,"lj/charmm/coul/long");            break;
    case LAMMPS_PAIR_CHCM: sprintf(s,"lj/charmm/coul/msm");             break;
    case LAMMPS_PAIR_CLS2: sprintf(s,"lj/class2");                      break;
    case LAMMPS_PAIR_C2CC: sprintf(s,"lj/class2/coul/cut");             break;
    case LAMMPS_PAIR_C2LC: sprintf(s,"lj/class2/coul/long");            break;
    case LAMMPS_PAIR_LJCB: sprintf(s,"lj/cubic");                       break;
    case LAMMPS_PAIR_LJCT: sprintf(s,"lj/cut");                         break;
    case LAMMPS_PAIR_LJCC: sprintf(s,"lj/cut/coul/cut");                break;
    case LAMMPS_PAIR_LJCD: sprintf(s,"lj/cut/coul/debye");              break;
    case LAMMPS_PAIR_LJCS: sprintf(s,"lj/cut/coul/dsf");                break;
    case LAMMPS_PAIR_LJCL: sprintf(s,"lj/cut/coul/long");               break;
    case LAMMPS_PAIR_LJCR: sprintf(s,"lj/cut/coul/long/cs");            break;
    case LAMMPS_PAIR_LJCM: sprintf(s,"lj/cut/coul/msm");                break;
    case LAMMPS_PAIR_LJDC: sprintf(s,"lj/cut/dipole/cut");              break;
    case LAMMPS_PAIR_LJDL: sprintf(s,"lj/cut/dipole/long");             break;
    case LAMMPS_PAIR_LJT4: sprintf(s,"lj/cut/tip4p/cut");               break;
    case LAMMPS_PAIR_LJL4: sprintf(s,"lj/cut/tip4p/long");              break;
    case LAMMPS_PAIR_LJVS: sprintf(s,"lj/expand");                      break;
    case LAMMPS_PAIR_GXLJ: sprintf(s,"lj/gromacs");                     break;
    case LAMMPS_PAIR_GXCL: sprintf(s,"lj/gromacs/coul/gromacs");        break;
    case LAMMPS_PAIR_LLCL: sprintf(s,"lj/long/coul/long");              break;
    case LAMMPS_PAIR_LLDL: sprintf(s,"lj/long/dipole/long");            break;
    case LAMMPS_PAIR_LLT4: sprintf(s,"lj/long/tip4p/long");             break;
    case LAMMPS_PAIR_SMLJ: sprintf(s,"lj/smooth");                      break;
    case LAMMPS_PAIR_SMLN: sprintf(s,"lj/smooth/linear");               break;
    case LAMMPS_PAIR_LJ96: sprintf(s,"lj96/cut");                       break;
    case LAMMPS_PAIR_LBPT: sprintf(s,"lubricate");                      break;
    case LAMMPS_PAIR_LBPL: sprintf(s,"lubricate/poly");                 break;
    case LAMMPS_PAIR_LBFD: sprintf(s,"lubricateU");                     break;
    case LAMMPS_PAIR_LBFP: sprintf(s,"lubricateU/poly");                break;
    case LAMMPS_PAIR_MEAM: sprintf(s,"meam");                           break;
    case LAMMPS_PAIR_MIEP: sprintf(s,"mie/cut");                        break;
    case LAMMPS_PAIR_MORS: sprintf(s,"morse");                          break;
    case LAMMPS_PAIR_NB3B: sprintf(s,"nb3b/harmonic");                  break;
    case LAMMPS_PAIR_NMPT: sprintf(s,"nm/cut");                         break;
    case LAMMPS_PAIR_NMCC: sprintf(s,"nm/cut/coul/cut");                break;
    case LAMMPS_PAIR_NMCL: sprintf(s,"nm/cut/coul/long");               break;
    case LAMMPS_PAIR_PREP: sprintf(s,"peri/eps");                       break;
    case LAMMPS_PAIR_PRLP: sprintf(s,"peri/lps");                       break;
    case LAMMPS_PAIR_PRPM: sprintf(s,"peri/pmb");                       break;
    case LAMMPS_PAIR_PRVS: sprintf(s,"peri/ves");                       break;
    case LAMMPS_PAIR_PLMR: sprintf(s,"polymorphic");                    break;
    case LAMMPS_PAIR_REAX: sprintf(s,"reax");                           break;
    case LAMMPS_PAIR_REBO: sprintf(s,"rebo");                           break;
    case LAMMPS_PAIR_RESQ: sprintf(s,"resquared");                      break;
    case LAMMPS_PAIR_SNAP: sprintf(s,"snap");                           break;
    case LAMMPS_PAIR_SOFT: sprintf(s,"soft");                           break;
    case LAMMPS_PAIR_SWPT: sprintf(s,"sw");                             break;
    case LAMMPS_PAIR_TABL: sprintf(s,"table");                          break;
    case LAMMPS_PAIR_TERS: sprintf(s,"tersoff");                        break;
    case LAMMPS_PAIR_TERM: sprintf(s,"tersoff/mod");                    break;
    case LAMMPS_PAIR_TERZ: sprintf(s,"tersoff/zbi");                    break;
    case LAMMPS_PAIR_TP4C: sprintf(s,"tip4p/cut");                      break;
    case LAMMPS_PAIR_TP4L: sprintf(s,"tip4p/long");                     break;
    case LAMMPS_PAIR_TRIL: sprintf(s,"tri/lj");                         break;
    case LAMMPS_PAIR_VASH: sprintf(s,"vashishta");                      break;
    case LAMMPS_PAIR_YUKW: sprintf(s,"yukawa");                         break;
    case LAMMPS_PAIR_YUKC: sprintf(s,"yukawa/colloid");                 break;
    case LAMMPS_PAIR_ZBLP: sprintf(s,"zbl");                            break;
    default:
      msg_error_f("unknown pair potential type (%d)",1,t);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
