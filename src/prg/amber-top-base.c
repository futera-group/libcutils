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
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* allocate new topology struct */
struct amber_top* amber_top_new(void) {
  struct amber_top *t=NULL;
  /* memory allocation */
  t = (struct amber_top*)malloc(sizeof(struct amber_top));
  if (!t)
    msg_error("cannot allocate memory for topology struct",1);
  /* initialization */
  t->version = NULL;
  t->title = NULL;
  t->pointers = NULL;
  t->atom_names = NULL;
  t->charge = NULL;
  t->atom_num = NULL;
  t->mass = NULL;
  t->atom_type_id = NULL;
  t->n_excl_atoms = NULL;
  t->nbond_prm_id = NULL;
  t->res_names = NULL;
  t->res_pointer = NULL;
  t->bond_frc_const = NULL;
  t->bond_eq_val = NULL;
  t->angle_frc_const = NULL;
  t->angle_eq_val = NULL;
  t->dihed_frc_const = NULL;
  t->dihed_period = NULL;
  t->dihed_phase = NULL;
  t->scale_ee = NULL;
  t->scale_nb = NULL;
  t->solty = NULL;
  t->lj_acoeff = NULL;
  t->lj_bcoeff = NULL;
  t->bond_h = NULL;
  t->bond_ah = NULL;
  t->angle_h = NULL;
  t->angle_ah = NULL;
  t->dihed_h = NULL;
  t->dihed_ah = NULL;
  t->excl_atom_list = NULL;
  t->hbond_acoeff = NULL;
  t->hbond_bcoeff = NULL;
  t->hbond_cut = NULL;
  t->atom_types = NULL;
  t->tree_chain_class = NULL;
  t->join_array = NULL;
  t->irotat = NULL;
  t->radii = NULL;
  t->screen = NULL;
  t->polar = NULL;
  t->box_pointer = NULL;
  t->box_natoms = NULL;
  t->box_params = NULL;
  t->radius_set = NULL;
  t->pol = 0;
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for topology struct

   t - pointer to topology struct */
void amber_top_free(struct amber_top *t) {
  if (t) {
    str_free(t->version);
    str_free(t->title);
    vec_sfree(t->atom_names,t->pointers[AMBER_POINTER_NATOM]);
    vec_sfree(t->res_names,t->pointers[AMBER_POINTER_NRES]);
    vec_sfree(t->atom_types,t->pointers[AMBER_POINTER_NATOM]);
    vec_sfree(t->tree_chain_class,t->pointers[AMBER_POINTER_NATOM]);
    vec_ifree(t->pointers);
    vec_ffree(t->charge);
    vec_ifree(t->atom_num);
    vec_ffree(t->mass);
    vec_ifree(t->atom_type_id);
    vec_ifree(t->n_excl_atoms);
    vec_ifree(t->nbond_prm_id);
    vec_ifree(t->res_pointer);
    vec_ffree(t->bond_frc_const);
    vec_ffree(t->bond_eq_val);
    vec_ffree(t->angle_frc_const);
    vec_ffree(t->angle_eq_val);
    vec_ffree(t->dihed_frc_const);
    vec_ffree(t->dihed_period);
    vec_ffree(t->dihed_phase);
    vec_ffree(t->scale_ee);
    vec_ffree(t->scale_nb);
    vec_ffree(t->solty);
    vec_ffree(t->lj_acoeff);
    vec_ffree(t->lj_bcoeff);
    vec_ifree(t->bond_h);
    vec_ifree(t->bond_ah);
    vec_ifree(t->angle_h);
    vec_ifree(t->angle_ah);
    vec_ifree(t->dihed_h);
    vec_ifree(t->dihed_ah);
    vec_ifree(t->excl_atom_list);
    vec_ffree(t->hbond_acoeff);
    vec_ffree(t->hbond_bcoeff);
    vec_ffree(t->hbond_cut);
    vec_ifree(t->join_array);
    vec_ifree(t->irotat);
    vec_ffree(t->radii);
    vec_ffree(t->screen);
    vec_ffree(t->polar);
    vec_ifree(t->box_pointer);
    vec_ifree(t->box_natoms);
    vec_ffree(t->box_params);
    str_free(t->radius_set);
    free(t);
    }
  }

/* -------------------------------------------------------------------------- */
