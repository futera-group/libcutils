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

#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* copy one topology struct to another

   t - pointer to topology struct from where data are read */
struct amber_top *amber_top_copy(struct amber_top *t) {
  struct amber_top *tn;
  struct amber_pointer p;
  /* data pointers */
  amber_pointer_init(&p,t->pointers);
  /* new topology structure */
  tn = amber_top_new();
  /* header */
  tn->version = str_copy_new(t->version);
  tn->title = str_copy_new(t->title);
  /* data pointers */
  tn->pointers = vec_icopy_new(t->pointers,AMBER_TOP_NPOINTER);
  /* atoms */
  tn->atom_names = vec_scopy_new(t->atom_names,p.natom,5);
  tn->charge = vec_fcopy_new(t->charge,p.natom);
  tn->mass = vec_fcopy_new(t->mass,p.natom);
  tn->atom_type_id = vec_icopy_new(t->atom_type_id,p.natom);
  tn->n_excl_atoms = vec_icopy_new(t->n_excl_atoms,p.natom);
  /* LJ parameteres */
  tn->nbond_prm_id = vec_icopy_new(t->nbond_prm_id,p.ntype*p.ntype);
  /* residui */
  tn->res_names = vec_scopy_new(t->res_names,p.nres,5);
  tn->res_pointer = vec_icopy_new(t->res_pointer,p.nres);
  /* bond force constants */
  tn->bond_frc_const = vec_fcopy_new(t->bond_frc_const,p.nbond);
  tn->bond_eq_val = vec_fcopy_new(t->bond_eq_val,p.nbond);
  /* angle force constants */
  tn->angle_frc_const = vec_fcopy_new(t->angle_frc_const,p.nangle);
  tn->angle_eq_val = vec_fcopy_new(t->angle_eq_val,p.nangle);
  /* dihedral force constants */
  tn->dihed_frc_const = vec_fcopy_new(t->dihed_frc_const,p.ndihed);
  tn->dihed_period = vec_fcopy_new(t->dihed_period,p.ndihed);
  tn->dihed_phase = vec_fcopy_new(t->dihed_phase,p.ndihed);
  tn->solty = vec_fcopy_new(t->solty,p.natyp);
  /* JL parameteres */
  tn->lj_acoeff = vec_fcopy_new(t->lj_acoeff,p.ntype*(p.ntype+1)/2);
  tn->lj_bcoeff = vec_fcopy_new(t->lj_bcoeff,p.ntype*(p.ntype+1)/2);
  /* bonds with / without hydrogens */
  tn->bond_h = vec_icopy_new(t->bond_h,3*p.nbondh);
  tn->bond_ah = vec_icopy_new(t->bond_ah,3*p.nbonda);
  /* angles with / without hydrogens */
  tn->angle_h = vec_icopy_new(t->angle_h,4*p.nangleh);
  tn->angle_ah = vec_icopy_new(t->angle_ah,4*p.nanglea);
  /* dihedrals with / without hydrogens */
  tn->dihed_h = vec_icopy_new(t->dihed_h,5*p.ndihedh);
  tn->dihed_ah = vec_icopy_new(t->dihed_ah,5*p.ndiheda);
  /* exluded atoms */
  tn->excl_atom_list = vec_icopy_new(t->excl_atom_list,p.next);
  /* hydrogen bonds */
  tn->hbond_acoeff = vec_fcopy_new(t->hbond_acoeff,p.nhbond);
  tn->hbond_bcoeff = vec_fcopy_new(t->hbond_bcoeff,p.nhbond);
  tn->hbond_cut = vec_fcopy_new(t->hbond_cut,p.nhbond);
  /* graph connectivities */
  tn->atom_types = vec_scopy_new(t->atom_types,p.natom,5);
  tn->tree_chain_class = vec_scopy_new(t->tree_chain_class,p.natom,5);
  tn->join_array = vec_icopy_new(t->join_array,p.natom);
  /* rotations */
  tn->irotat = vec_icopy_new(t->irotat,p.natom);
  /* atom radii and screening */
  tn->radii = vec_fcopy_new(t->radii,p.natom);
  tn->screen = vec_fcopy_new(t->screen,p.natom);
  /* polarizabilities */
  tn->polar = (t->polar ? vec_fcopy_new(t->polar,p.natom) : NULL);
  /* solvent */
  tn->box_pointer = vec_icopy_new(t->box_pointer,3);
  tn->box_natoms = vec_icopy_new(t->box_natoms,t->box_pointer[AMBER_BOX_NMOL]);
  tn->box_params = vec_fcopy_new(t->box_params,4);
  return(tn);
  }

/* -------------------------------------------------------------------------- */
