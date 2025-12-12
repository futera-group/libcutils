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
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* get angular vector of specific cartesian gaussian type orbital
 
   t  - type of the shell (s,p,d,...)
   s  - shell basis function ID
   id - angular vector (output) */
void gto_ang_vec(short t, short s, unsigned *id) {
  /* divide SP shell */
  if (t==BASIS_SHELL_SP) {
    t = (s>0 ? BASIS_SHELL_P : BASIS_SHELL_S);
    s = (s>0 ? s-1 : s);
    }
  /* set cartesian ID */
  switch (t) {
    case BASIS_SHELL_S:
      id[0] = 0;
      id[1] = 0;
      id[2] = 0;
      break;
    case BASIS_SHELL_P:
      switch (s) {
        case 0: 
          id[0] = 1;
          id[1] = 0;
          id[2] = 0;
          break;
        case 1:
          id[0] = 0;
          id[1] = 1;
          id[2] = 0;
          break;
        case 2:
          id[0] = 0;
          id[1] = 0;
          id[2] = 1;
          break;
        }
      break;
    case BASIS_SHELL_Dc:
      switch (s) {
        case 0:
          id[0] = 2;
          id[1] = 0;
          id[2] = 0;
          break;
        case 1:
          id[0] = 0;
          id[1] = 2;
          id[2] = 0;
          break;
        case 2:
          id[0] = 0;
          id[1] = 0;
          id[2] = 2;
          break;
        case 3:
          id[0] = 1;
          id[1] = 1;
          id[2] = 0;
          break;
        case 4:
          id[0] = 1;
          id[1] = 0;
          id[2] = 1;
          break;
        case 5:
          id[0] = 0;
          id[1] = 1;
          id[2] = 1;
          break;
        }
      break;
    default:
      msg_error_f("unsupported GTO type: %d/%d",1,t,s);
    }
  }

/* -------------------------------------------------------------------------- */

/* get name of specified angular momemntum

   a - the angular momentum */
char gto_ang_name(unsigned a) {
  char lab = 'X';
  switch (a) {
    case 0: lab = 'S'; break;
    case 1: lab = 'P'; break;
    case 2: lab = 'D'; break;
    case 3: lab = 'F'; break;
    }
  return(lab);
  }

/* -------------------------------------------------------------------------- */
