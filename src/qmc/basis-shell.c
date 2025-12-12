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
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* return number of shells in the basis set
 
   b - pointer to basis set data struct */
unsigned basis_shell_num(struct basis *b) {
  unsigned i,ns = 0;
  for (i=0; i<b->n_centers; i++)
    ns += b->center[i].n_shells;
  return(ns);
  }

/* return number of basis functions in the shell
 
   s - shell ID */
unsigned basis_shell_nbfce(short s) {
  switch (s) {
    case BASIS_SHELL_S:  return(1);
    case BASIS_SHELL_SP: return(4);
    case BASIS_SHELL_P:  return(3);
    case BASIS_SHELL_Dp: return(5);
    case BASIS_SHELL_Dc: return(6);
    case BASIS_SHELL_Fp: return(7);
    case BASIS_SHELL_Fc: return(10);
    case BASIS_SHELL_Gp: return(9);
    case BASIS_SHELL_Gc: return(15);
    case BASIS_SHELL_Hp: return(11);
    case BASIS_SHELL_Hc: return(21);
    case BASIS_SHELL_Ip: return(13);
    case BASIS_SHELL_Ic: return(28);
    case BASIS_SHELL_Jp: return(15);
    case BASIS_SHELL_Jc: return(36);
    }
  return(0);
  }

/* -------------------------------------------------------------------------- */

/* convert shell id to string name 

   id   - id of the shell */
char* basis_shell_name(short id) {
  static char name[5]="\0";
  switch (id) {
    case BASIS_SHELL_X:  sprintf(name,"X");  break;
    case BASIS_SHELL_S:  sprintf(name,"S");  break;
    case BASIS_SHELL_SP: sprintf(name,"SP"); break;
    case BASIS_SHELL_P:  sprintf(name,"P");  break;
    case BASIS_SHELL_Dp: sprintf(name,"D");  break;
    case BASIS_SHELL_Dc: sprintf(name,"D");  break;
    case BASIS_SHELL_Fp: sprintf(name,"F");  break;
    case BASIS_SHELL_Fc: sprintf(name,"F");  break;
    case BASIS_SHELL_Gp: sprintf(name,"G");  break;
    case BASIS_SHELL_Gc: sprintf(name,"G");  break;
    case BASIS_SHELL_Hp: sprintf(name,"H");  break;
    case BASIS_SHELL_Hc: sprintf(name,"H");  break;
    case BASIS_SHELL_Ip: sprintf(name,"I");  break;
    case BASIS_SHELL_Ic: sprintf(name,"I");  break;
    case BASIS_SHELL_Jp: sprintf(name,"J");  break;
    case BASIS_SHELL_Jc: sprintf(name,"J");  break;
    default:             sprintf(name,"X");  break;
    }
  return(name);
  }

/* convert shell id to string function name 

   t - shell type
   n - id of the function */
char* basis_shell_name_fce(short t, unsigned n) {
  static char name[5]="\0";
  char label[25]="\0";
  /* label */
  if (t==BASIS_SHELL_SP) {
    if (n)
      sprintf(label,"%s",basis_shell_name(BASIS_SHELL_P));
    else
      sprintf(label,"%s",basis_shell_name(BASIS_SHELL_S));
    }
  else
    sprintf(label,"%s",basis_shell_name(t));
  /* number */
  if (t==BASIS_SHELL_S || (t==BASIS_SHELL_SP && !n))
    sprintf(name,"%s",label);
  else if (t==BASIS_SHELL_SP)
    sprintf(name,"%s%d",label,n);
  else
    sprintf(name,"%s%d",label,n+1);
  return(name);
  }

/* -------------------------------------------------------------------------- */

/* return ID of the shell which name is given

   name - name of the shell
   pure - pure of cartesian functions */
short basis_shell_id(char *name, short pure) {
  unsigned id = BASIS_SHELL_X;
  if (str_compare(name,"S"))
    id = BASIS_SHELL_S;
  else if (str_compare(name,"SP"))
    id = BASIS_SHELL_SP;
  else if (str_compare(name,"P"))
    id = BASIS_SHELL_P;
  else if (str_compare(name,"D"))
    id = (pure ? BASIS_SHELL_Dp : BASIS_SHELL_Dc);
  else if (str_compare(name,"F"))
    id = (pure ? BASIS_SHELL_Fp : BASIS_SHELL_Fc);
  else if (str_compare(name,"G"))
    id = (pure ? BASIS_SHELL_Gp : BASIS_SHELL_Gc);
  else if (str_compare(name,"H"))
    id = (pure ? BASIS_SHELL_Hp : BASIS_SHELL_Hc);
    else if (str_compare(name,"I"))
    id = (pure ? BASIS_SHELL_Ip : BASIS_SHELL_Ic);
  else if (str_compare(name,"J"))
    id = (pure ? BASIS_SHELL_Jp : BASIS_SHELL_Jc);
  return(id);
  }

/* -------------------------------------------------------------------------- */

/* return angular momentum of specified shell 
 
   s - shell ID */
unsigned basis_shell_ang_mom(short s) {
  switch (s) {
    case BASIS_SHELL_S:  return(0);
    case BASIS_SHELL_SP: return(1);
    case BASIS_SHELL_P:  return(1);
    case BASIS_SHELL_Dp: return(2);
    case BASIS_SHELL_Dc: return(2);
    case BASIS_SHELL_Fp: return(3);
    case BASIS_SHELL_Fc: return(3);
    case BASIS_SHELL_Gp: return(4);
    case BASIS_SHELL_Gc: return(4);
    case BASIS_SHELL_Hp: return(5);
    case BASIS_SHELL_Hc: return(5);
    case BASIS_SHELL_Ip: return(6);
    case BASIS_SHELL_Ic: return(6);
    case BASIS_SHELL_Jp: return(7);
    case BASIS_SHELL_Jc: return(7);
    }
  return(0);
  }

/* -------------------------------------------------------------------------- */
