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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of run type 
 
   key - run type keyword */
short cp2k_run_type_id(char *key) {
  char s[80] = "\n";
  str_lowcase_copy(key,s);
  if (str_compare(s,"band"))
    return(RUN_TYPE_BAND);
  if (str_compare(s,"bsse"))
    return(RUN_TYPE_BSSE);
  if (str_compare(s,"cell_opt"))
    return(RUN_TYPE_COPT);
  if (str_compare(s,"debug"))
    return(RUN_TYPE_DEBG);
  if (str_compare(s,"driver"))
    return(RUN_TYPE_DRVR);
  if (str_compare(s,"ehrenfest_dyn"))
    return(RUN_TYPE_EFST);
  if (str_compare(s,"electronic_spectra"))
    return(RUN_TYPE_SPCT);
  if (str_compare(s,"energy"))
    return(RUN_TYPE_WOPT);
  if (str_compare(s,"energy_force"))
    return(RUN_TYPE_EFRC);
  if (str_compare(s,"geometry_optimization"))
    return(RUN_TYPE_GOPT);
  if (str_compare(s,"geo_opt"))
    return(RUN_TYPE_GOPT);
  if (str_compare(s,"linear_response"))
    return(RUN_TYPE_LINR);
  if (str_compare(s,"lr"))
    return(RUN_TYPE_LINR);
  if (str_compare(s,"mc"))
    return(RUN_TYPE_MNCR);
  if (str_compare(s,"md"))
    return(RUN_TYPE_MDYN);
  if (str_compare(s,"molecular_dynamics"))
    return(RUN_TYPE_MDYN);
  if (str_compare(s,"montecarlo"))
    return(RUN_TYPE_MNCR);
  if (str_compare(s,"none"))
    return(RUN_TYPE_NONE);
  if (str_compare(s,"normal_modes"))
    return(RUN_TYPE_NMOD);
  if (str_compare(s,"pint"))
    return(RUN_TYPE_PINT);
  if (str_compare(s,"rt_propagation"))
    return(RUN_TYPE_RTPR);
  if (str_compare(s,"spectra"))
    return(RUN_TYPE_SPCT);
  if (str_compare(s,"tamc"))
    return(RUN_TYPE_TAMC);
  if (str_compare(s,"tmc"))
    return(RUN_TYPE_TRMC);
  if (str_compare(s,"vibrational_analysis"))
    return(RUN_TYPE_NMOD);
  if (str_compare(s,"wavefunction_optmization"))
    return(RUN_TYPE_WOPT);
  if (str_compare(s,"wfn_opt"))
    return(RUN_TYPE_WOPT);
  msg_error_f("unknown run type keyword \"%s\"",1,key);
  return(0);
  }

/* convert internal run type ID to corresponding keyword
 
   key - storage for the keyword
   id  - the internal ID */
void cp2k_run_type_name(char *key, short id) {
  switch (id) {
    case RUN_TYPE_BAND: sprintf(key,"band");           break;
    case RUN_TYPE_BSSE: sprintf(key,"bsse");           break;
    case RUN_TYPE_COPT: sprintf(key,"cell_opt");       break;
    case RUN_TYPE_DEBG: sprintf(key,"debug");          break;
    case RUN_TYPE_DRVR: sprintf(key,"driver");         break;
    case RUN_TYPE_EFST: sprintf(key,"ehrenfest_dyn");  break;
    case RUN_TYPE_SPCT: sprintf(key,"spectra");        break;
    case RUN_TYPE_WOPT: sprintf(key,"energy");         break;
    case RUN_TYPE_EFRC: sprintf(key,"energy_force");   break;
    case RUN_TYPE_GOPT: sprintf(key,"geo_opt");        break;
    case RUN_TYPE_LINR: sprintf(key,"lr");             break;
    case RUN_TYPE_MNCR: sprintf(key,"mc");             break;
    case RUN_TYPE_MDYN: sprintf(key,"md");             break;
    case RUN_TYPE_NONE: sprintf(key,"none");           break;
    case RUN_TYPE_NMOD: sprintf(key,"normal_mode");    break;
    case RUN_TYPE_PINT: sprintf(key,"pint");           break;
    case RUN_TYPE_RTPR: sprintf(key,"rt_propagation"); break;
    case RUN_TYPE_TAMC: sprintf(key,"tamc");           break;
    case RUN_TYPE_TRMC: sprintf(key,"tmc");            break;
    default:
      msg_error_f("unknown internal run-type ID (%d)",1,id);
    }                      
  }

/* -------------------------------------------------------------------------- */
