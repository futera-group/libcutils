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

#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* return list of atoms bonded to specified atom
 
   c - pointer to acf molecular data struct
   a - ID of the atom
   n - number of bonds (output) */
unsigned* mol_acf_type_gaff_atom_bonds(struct acf_mol *c,
  unsigned a, unsigned *n) {
  unsigned ib,*b = NULL;
  struct queue *q;
  q = queue_alloc();
  /* find bonds */
  for (ib=0; ib<c->n_bonds; ib++) {
    if (c->bond[ib][0]==a)
      queue_uadd(q,c->bond[ib][1]);
    if (c->bond[ib][1]==a)
      queue_uadd(q,c->bond[ib][0]);
    }
  /* create atom ID array */
  (*n) = q->num;
  b = vec_ualloc(q->num);
  ib = 0;
  while (q->num) {
    queue_uget(q,&(b[ib]));
    ib++;
    }
  /* clean memory */
  queue_free(q);
  return(b);
  }

/* Try to assign Amber GAFF type to hydrogen atoms

   c - pointer to acf molecular data struct
   t - atomic type (output)
   b - array of bonded atoms 
   n - number of bonds */
void mol_acf_type_gaff_set_h(struct acf_mol *c, char *t, unsigned *b, unsigned n) {
  sprintf(t,"h");
  if (n==1) {
    switch (c->atom[b[0]].num) {
      case  6: sprintf(t,"hc"); break;
      case  8: sprintf(t,"ho"); break;
      case 15: sprintf(t,"hp"); break;
      case 16: sprintf(t,"hs"); break;
      }
    }
  }

/* Try to assign Amber GAFF type to carbon atoms

   c - pointer to acf molecular data struct
   t - atomic type (output)
   b - array of bonded atoms 
   n - number of bonds */
void mol_acf_type_gaff_set_c(struct acf_mol *c, char *t, unsigned *b, unsigned n) {
  sprintf(t,"c");
  switch (n) {
    case 2: sprintf(t,"c1"); break;
    case 3: sprintf(t,"c2"); break;
    case 4: sprintf(t,"c3"); break;
    }
  }

/* Try to assign Amber GAFF type to boron atoms

   c - pointer to acf molecular data struct
   t - atomic type (output)
   b - array of bonded atoms 
   n - number of bonds */
void mol_acf_type_gaff_set_b(struct acf_mol *c, char *t, unsigned *b, unsigned n) {
  sprintf(t,"b");
  }

/* Try to assign Amber GAFF type to nitrogen atoms

   c - pointer to acf molecular data struct
   t - atomic type (output)
   b - array of bonded atoms 
   n - number of bonds */
void mol_acf_type_gaff_set_n(struct acf_mol *c, char *t, unsigned *b, unsigned n) {
  unsigned i,nh=0,no=0;
  sprintf(t,"n");
  switch (n) {
    case 3:
      for (i=0; i<n; i++) {
        if (c->atom[b[i]].num==1)
          nh++;
        if (c->atom[b[i]].num==8)
          no++;
        }
      if (no==2)
        sprintf(t,"no");
      else if (nh==2)
        sprintf(t,"nh");
      else
        sprintf(t,"n3");
      break;
    case 4: 
      sprintf(t,"n4");
      break;
    }
  }

/* Try to assign Amber GAFF type to oxygen atoms

   c - pointer to acf molecular data struct
   t - atomic type (output)
   b - array of bonded atoms 
   n - number of bonds */
void mol_acf_type_gaff_set_o(struct acf_mol *c, char *t, unsigned *b, unsigned n) {
  sprintf(t,"o");
  switch (n) {
    case 1: 
      if (c->atom[b[0]].num==1)
        sprintf(t,"oh");
      break;
    case 2: 
      if (c->atom[b[0]].num==1 && c->atom[b[1]].num==1)
        sprintf(t,"ow");
      else if (c->atom[b[0]].num==1 || c->atom[b[1]].num==1)
        sprintf(t,"oh");
      else
        sprintf(t,"os");
      break;
    }
  }

/* Try to assign Amber GAFF types to all atoms

   c - pointer to acf molecular data struct */
void mol_acf_type_gaff(struct acf_mol *c) {
  unsigned ia,nb,*b;
  char type[5];
  for (ia=0; ia<c->n_atoms; ia++) {
    b = mol_acf_type_gaff_atom_bonds(c,ia,&nb);
    type[0] = '\0';
    switch (c->atom[ia].num) {
      case 1: mol_acf_type_gaff_set_h(c,type,b,nb); break;
      case 5: mol_acf_type_gaff_set_b(c,type,b,nb); break;
      case 6: mol_acf_type_gaff_set_c(c,type,b,nb); break;
      case 7: mol_acf_type_gaff_set_n(c,type,b,nb); break;
      case 8: mol_acf_type_gaff_set_o(c,type,b,nb); break;
      }
    if (type[0]) {
      c->atom[ia].type = str_free(c->atom[ia].type);
      c->atom[ia].type = str_copy_new(type);
      }
    vec_ufree(b);
    }
  }

/* -------------------------------------------------------------------------- */
