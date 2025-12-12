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

#ifndef ZF_LIB_MOL_MOLEC_H
#define ZF_LIB_MOL_MOLEC_H

/* symbolic constants */
#define TREE_MARK_M  1     /* main chain atom */
#define TREE_MARK_E  2     /* ending atom -1 connection */
#define TREE_MARK_S  3     /* side chain atom - 2 connections */
#define TREE_MARK_B  4     /* branch atom - 3 connections */
#define TREE_MARK_3  5     /* atom with 4 connections */
#define TREE_MARK_4  6     /* atom with 5 connections */
#define TREE_MARK_5  7     /* atom with 6 connections */
#define TREE_MARK_6  8     /* atom with 7 connections */

/* tree atomic data struct */
struct mol_tree_atom {
  unsigned id;             /* atom ID */
  short tree;              /* tree mark */
  };

/* allocate memory for new tree atomic data structure */
struct mol_tree_atom* mol_tree_atom_new(void);
/* free memory allocated for a tree atomic data structure */
void mol_tree_atom_free(struct mol_tree_atom *t);

/* convert bond molecular representation to graph tree */
struct tree *mol_tree_set(unsigned**, unsigned, unsigned, unsigned***,
  unsigned*);

/* General molecular format */
#include "mol-gen.h"

/* Amber AC format */
#include "mol-acf.h"
/* Amber Prep format */
#include "mol-apc.h"
/* CIF format */
#include "mol-cif.h"
/* Gaussian CUBE format */
#include "mol-cub.h"
/* PDB format */
#include "mol-pdb.h"
/* XYZ format */
#include "mol-xyz.h"
/* Z-Matrix format */
#include "mol-zmt.h"

#endif
