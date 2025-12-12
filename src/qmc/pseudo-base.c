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

#include <stdlib.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/vector.h>
#include "qmc/pseudo.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for pseudopotential data */
struct pseudo *pseudo_new(void) {
  struct pseudo *p = NULL;
  /* memory allocation */
  p = (struct pseudo*)malloc(sizeof(struct pseudo));
  if (!p)
    msg_error("cannot allocate memory for pseudopotential data",1);
  /* initialization */
  p->atom_num = 0;
  p->val_num = 0;
  p->xc_num = 0;
  p->xc_slater = 0;
  p->xc_lda_c = 0;
  p->xc_gc_e = 0;
  p->xc_gc_c = 0;
  p->pp_type = 0;
  p->time_stamp = NULL;
  p->state_core = 0;
  p->state_val = 0;
  p->mesh_num = 0;
  p->en_pot = 0.0;
  p->en_tot = 0.0;
  p->el_conf_n = NULL;
  p->el_conf_l = NULL;
  p->el_conf_occ = NULL;
  p->mt_num = 0;
  p->mt_pp_n = NULL;
  p->mt_pp_l = NULL;
  p->mt_pp_r = NULL;
  p->mt_pp_e = NULL;
  p->dat_pot = NULL;
  p->dat_wfce = NULL;
  p->dat_den = NULL;
  return(p);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for pseudopotential data

   p - pointer to pseudopotential data struct */
void pseudo_free(struct pseudo *p) {
  if (p) {
    p->el_conf_n = vec_ufree(p->el_conf_n);
    p->el_conf_l = vec_ufree(p->el_conf_l);
    p->el_conf_occ = vec_ffree(p->el_conf_occ);
    p->mt_pp_n = vec_ufree(p->mt_pp_n);
    p->mt_pp_l = vec_sifree(p->mt_pp_l);
    p->mt_pp_r = vec_ffree(p->mt_pp_r);
    p->mt_pp_e = vec_ffree(p->mt_pp_e);
    p->dat_pot = mat_ffree(p->dat_pot,p->mesh_num);
    p->dat_wfce = mat_ffree(p->dat_wfce,p->mesh_num);
    p->dat_den = mat_ffree(p->dat_den,p->mesh_num);
    free(p);
    }
  }

/* -------------------------------------------------------------------------- */
