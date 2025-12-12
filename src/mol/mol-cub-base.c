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
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate new cube molecule struct */
struct cub_mol *mol_cub_new(void) {
  /* memory allocation */
  struct cub_mol *c = NULL;
  c = (struct cub_mol*)malloc(sizeof(struct cub_mol));
  if (!c)
    msg_error("cannot allocate memory for cube molecular data struct",1);
  /* initialization */
  c->title = NULL;
  c->desc = NULL;
  c->n_atoms = 0;
  c->atom = NULL;
  c->grid = NULL;
  c->data = NULL;
  return(c);
  }

/* allocate array of cube atom data structures 

   n - number of atoms */
struct cub_atom *mol_cub_atom_new(unsigned n) {
  unsigned i,j;
  /* memory allocation */
  struct cub_atom *c = NULL;
  c = (struct cub_atom*)malloc(n*sizeof(struct cub_atom));
  if (!c)
    msg_error("cannot allocate memory for array of cube atomr data struct",1);
  /* initialization */
  for (i=0; i<n; i++) {
    c[i].num = 0;
    for (j=0; j<3; j++)
      c[i].coord[j] = 0.0;
    }
  return(c);
  }

/* allocate new cube grid parameters struct */
struct cub_grid *mol_cub_grid_new(void) {
  unsigned i,j;
  struct cub_grid *g = NULL;
  /* memory allocation */
  g = (struct cub_grid*)malloc(sizeof(struct cub_grid));
  if (!g)
    msg_error("cannot allocate memory for cube grid parameters struct",1);
  /* initialization */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++)
      g->vector[i][j] = 0.0;
    g->origin[i] = 0.0;
    g->range[i] = 0;
    }
  return(g);
  }

/* allocate new grid data array
 
   g - pointer to grid parameters struct */
double*** mol_cub_data_new(struct cub_grid *g) {
  unsigned i,j,k;
  double ***d = NULL;
  d = (double***)malloc(g->range[0]*sizeof(double**));
  if (!d)
    msg_error("cannot allocate new cube data array",1);
  for (i=0; i<g->range[0]; i++) {
    d[i] = (double**)malloc(g->range[1]*sizeof(double*));
    if (!d[i])
      msg_error("cannot allocate new cube data array",1);
    for (j=0; j<g->range[1]; j++) {
      d[i][j] = (double*)malloc(g->range[2]*sizeof(double));
      if (!d[i][j])
        msg_error("cannot allocate new cube data array",1);
      for (k=0; k<g->range[2]; k++)
        d[i][j][k] = 0.0;
      }
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for cube molecule struct

   c - pointer to xyz molecular data struct */
void mol_cub_free(struct cub_mol *c) {
  if (c) {
    str_free(c->title);
    str_free(c->desc);
    mol_cub_atom_free(c->atom);
    mol_cub_grid_free(c->grid);
    mol_cub_data_free(c->grid,c->data);
    free(c);
    }
  }

/* free memory allocated for array of atom data structures 

   c - array of atoms */
void mol_cub_atom_free(struct cub_atom *c) {
  if (c) 
    free(c);
  }

/* free memory allocated for grid parameter data structure

   c - grid parameter data */
void mol_cub_grid_free(struct cub_grid *c) {
  if (c) 
    free(c);
  }

/* free memory allocated for grid data array

   g - grid parameters
   d - grid data array */
void mol_cub_data_free(struct cub_grid *g, double ***d) {
  unsigned i,j;
  if (d) {
    for (i=0; i<g->range[0]; i++)
      if (d[i]) {
        for (j=0; j<g->range[1]; j++) {
          if (d[i][j])
            free(d[i][j]);
          }
        free(d[i]);
        }
    free(d);
    }
  }

/* -------------------------------------------------------------------------- */

/* create copy of cube grid parameters struct
 
   c - pointer to the cube grid data struct */
struct cub_grid *mol_cub_grid_copy_new(struct cub_grid *c) {
  unsigned i,j;
  struct cub_grid *g = NULL;
  g = mol_cub_grid_new();
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++)
      g->vector[i][j] = c->vector[i][j];
    g->origin[i] = c->origin[i];
    g->range[i] = c->range[i];
    }
  return(g);
  }

/* -------------------------------------------------------------------------- */
