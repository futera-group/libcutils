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

#include <math.h>
#include <stdlib.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "prg/dlpoly.h"
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* check if the provided potential is unique and add to queue
 
   t - the potential data
   q - the queue with the potentials */
unsigned lammps_conv_pot_find(struct lammps_pot *t, struct queue *q) {
  unsigned i,t_id;
  short found;
  struct lammps_pot *r;
  struct qdata *p;
  for (p=q->data, t_id=0, found=0; p; p=p->q_next, t_id++) {
    r = p->q_data;
    if (r->type==t->type && r->n_coeffs==t->n_coeffs) {
      for (i=0; i<r->n_coeffs; i++)
        if (fabs(r->coeff[i]-t->coeff[i])>1.0e-5) 
          break;
      if (i>=r->n_coeffs) {
        found = 1;
        break;
        }
      }
    }
  if (!found)
    queue_add(q,t);
  else {
    t->coeff = vec_ffree(t->coeff);
    free(t);
    }
  return(t_id);
  }

/* convert potential definitions saved in queue to array in LAMMPS data storage
 
   p - the potential array
   n - number of the potential definitions (output)
   q - the queue with the potentials */
void lammps_conv_pot_save(struct lammps_pot **p, unsigned *n, struct queue *q) {
  unsigned i;
  struct lammps_pot *t;
  (*n) = q->num;
  (*p) = lammps_pot_new(*n);
  for (i=0; i<(*n); i++) {
    t = queue_get(q);
    ((*p)[i]).type = t->type;
    ((*p)[i]).style = t->style;
    ((*p)[i]).n_ids = t->n_ids;
    ((*p)[i]).id = vec_ucopy_new(t->id,t->n_ids);
    t->id = vec_ufree(t->id);
    ((*p)[i]).n_coeffs = t->n_coeffs;
    ((*p)[i]).coeff = vec_fcopy_new(t->coeff,t->n_coeffs);
    t->coeff = vec_ffree(t->coeff);
    free(t);
    }
  }

/* convert pair potentials from DL_POLY format to LAMMPS
 
   d - LAMMPS structure & force-field data
   f - DL_POLY force-field data */
void lammps_from_dlpoly_pairs(struct lammps_dat *d, struct dlpoly_fld *f) {
  unsigned i,j;
  /* check presence of VDW parameters */
  if (f->n_vdws && f->vdw) {
    d->pot_pair = lammps_pot_new(d->n_atom_types);
    for (i=0; i<d->n_atom_types; i++) {
      d->pot_pair[i].style = LAMMPS_PAIR_LJCT;
      d->pot_pair[i].n_coeffs = 2;
      d->pot_pair[i].coeff = vec_falloc(d->pot_pair[i].n_ids);
      d->pot_pair[i].coeff[0] = 0.0;
      d->pot_pair[i].coeff[1] = 0.0;
      /* epsilon / sigma coefficient */
      for (j=0; j<f->n_vdws; j++)
        if (str_compare(f->vdw[j].at1,d->type[i].name) &&
            str_compare(f->vdw[j].at2,d->type[i].name)) {
          switch (f->vdw[j].pot) {
            case DLPOLY_PAIR_12: 
              if (fabs(f->vdw[j].v1)>1.0e-5 && fabs(f->vdw[j].v2)>1.0e-5) {
                d->pot_pair[i].coeff[0] = 
                  (f->vdw[j].v2*f->vdw[j].v2)/(4.0*f->vdw[j].v1);
                d->pot_pair[i].coeff[1] = 
                  pow(f->vdw[j].v1/f->vdw[j].v2,1.0/6.0);
                }
              break;
            case DLPOLY_PAIR_LJ: 
              d->pot_pair[i].coeff[0] = f->vdw[j].v1;
              d->pot_pair[i].coeff[1] = f->vdw[j].v2;
              break;
            }
          break;
          }
      }
    }
  }

/* convert dihedral-angle definitions from DL_POLY format to LAMMPS
 
   d - LAMMPS structure & force-field data
   f - DL_POLY force-field data */
void lammps_from_dlpoly_dihedrals(struct lammps_dat *d, struct dlpoly_fld *f) {
  unsigned i,j,k,id;
  struct lammps_pot *t,*a;
  struct queue *qt,*qa;
  /* initialization */
  qt = queue_alloc();
  qa = queue_alloc();
  id = 0;
  /* search all defined bonds */
  for (i=0; i<f->n_molecules; i++)
    for (j=0; j<f->mol[i].n_molecules; j++) {
      /* dihedral angles */
      for (k=0; k<f->mol[i].n_dihedrals; k++) {
        /* atom IDs */
        a = lammps_pot_new(1);
        a->n_ids = 4;
        a->id = vec_ualloc(a->n_ids);
        a->id[0] = f->mol[i].dihed[k].id1 + id;
        a->id[1] = f->mol[i].dihed[k].id2 + id;
        a->id[2] = f->mol[i].dihed[k].id3 + id;
        a->id[3] = f->mol[i].dihed[k].id4 + id;
        /* coefficients */
        t = lammps_pot_new(1);
        switch (f->mol[i].dihed[k].pot) {
          case DLPOLY_DIHE_CS: /* cosine potential */
            t->style = LAMMPS_DIHE_CHRM;
            t->n_coeffs = 4;
            t->coeff = vec_falloc(t->n_coeffs);
            t->coeff[0] = f->mol[i].dihed[k].v1;
            t->coeff[1] = f->mol[i].dihed[k].v3;
            t->coeff[2] = f->mol[i].dihed[k].v2;
            t->coeff[3] = 0.0;
            break;
          default:
            msg_warn_f("dihedral potential \"%s\" not converted",
              dlpoly_prm_dihed_pot_sym(f->mol[i].dihed[k].pot));
          }
        /* potential ID */
        a->type = lammps_conv_pot_find(t,qt);
        queue_add(qa,a);
        }
      /* atom-ID offset */
      for (k=0; k<f->mol[i].n_atoms; k++)
        id += f->mol[i].atom[k].n_atoms;
      }
  /* convert potentials */
  lammps_conv_pot_save(&(d->dihed),&(d->n_dihedrals),qa);
  lammps_conv_pot_save(&(d->pot_dihed),&(d->n_dihed_types),qt);
  /* clean memory */
  queue_free(qt);
  queue_free(qa);
  }

/* convert valence-angle definitions from DL_POLY format to LAMMPS
 
   d - LAMMPS structure & force-field data
   f - DL_POLY force-field data */
void lammps_from_dlpoly_angles(struct lammps_dat *d, struct dlpoly_fld *f) {
  unsigned i,j,k,id;
  struct lammps_pot *t,*a;
  struct queue *qt,*qa;
  /* initialization */
  qt = queue_alloc();
  qa = queue_alloc();
  id = 0;
  /* search all defined bonds */
  for (i=0; i<f->n_molecules; i++)
    for (j=0; j<f->mol[i].n_molecules; j++) {
      /* valence angles */
      for (k=0; k<f->mol[i].n_angles; k++) {
        /* atom IDs */
        a = lammps_pot_new(1);
        a->n_ids = 3;
        a->id = vec_ualloc(a->n_ids);
        a->id[0] = f->mol[i].angle[k].id1 + id;
        a->id[1] = f->mol[i].angle[k].id2 + id;
        a->id[2] = f->mol[i].angle[k].id3 + id;
        /* coefficients */
        t = lammps_pot_new(1);
        switch (f->mol[i].angle[k].pot) {
          case DLPOLY_BOND_HM: /* harmonic potential */
            t->style = LAMMPS_ANGL_HARM;
            t->n_coeffs = 2;
            t->coeff = vec_falloc(t->n_coeffs);
            t->coeff[0] = f->mol[i].angle[k].v1/2.0;
            t->coeff[1] = f->mol[i].angle[k].v2;
            break;
          default:
            msg_warn_f("angle potential \"%s\" not converted",
              dlpoly_prm_angle_pot_sym(f->mol[i].angle[k].pot,0));
          }
        /* potential ID */
        a->type = lammps_conv_pot_find(t,qt);
        queue_add(qa,a);
        }
      /* atom-ID offset */
      for (k=0; k<f->mol[i].n_atoms; k++)
        id += f->mol[i].atom[k].n_atoms;
      }
  /* convert potentials */
  lammps_conv_pot_save(&(d->angle),&(d->n_angles),qa);
  lammps_conv_pot_save(&(d->pot_angle),&(d->n_angle_types),qt);
  /* clean memory */
  queue_free(qt);
  queue_free(qa);
  }

/* convert bond definitions from DL_POLY format to LAMMPS
 
   d - LAMMPS structure & force-field data
   f - DL_POLY force-field data */
void lammps_from_dlpoly_bonds(struct lammps_dat *d, struct dlpoly_fld *f) {
  unsigned i,j,k,id;
  struct lammps_pot *t,*b;
  struct queue *qt,*qb;
  /* initialization */
  qt = queue_alloc();
  qb = queue_alloc();
  id = 0;
  /* search all defined bonds */
  for (i=0; i<f->n_molecules; i++)
    for (j=0; j<f->mol[i].n_molecules; j++) {
      /* interatomic bonds */
      for (k=0; k<f->mol[i].n_bonds; k++) {
        /* atom IDs */
        b = lammps_pot_new(1);
        b->n_ids = 2;
        b->id = vec_ualloc(b->n_ids);
        b->id[0] = f->mol[i].bond[k].id1 + id;
        b->id[1] = f->mol[i].bond[k].id2 + id;
        /* coefficients */
        t = lammps_pot_new(1);
        switch (f->mol[i].bond[k].pot) {
          case DLPOLY_BOND_HM: /* harmonic potential */
            t->style = LAMMPS_BOND_HARM;
            t->n_coeffs = 2;
            t->coeff = vec_falloc(t->n_coeffs);
            t->coeff[0] = f->mol[i].bond[k].v1/2.0;
            t->coeff[1] = f->mol[i].bond[k].v2;
            break;
          case DLPOLY_BOND_MR: /* Morse potential */
            t->style = LAMMPS_BOND_MORS;
            t->n_coeffs = 3;
            t->coeff = vec_falloc(t->n_coeffs);
            t->coeff[0] = f->mol[i].bond[k].v1;
            t->coeff[1] = f->mol[i].bond[k].v3;
            t->coeff[2] = f->mol[i].bond[k].v2;
            break;
          case DLPOLY_BOND_QR: /* quartic potential */
            t->style = LAMMPS_BOND_CLS2;
            t->n_coeffs = 4;
            t->coeff = vec_falloc(t->n_coeffs);
            t->coeff[0] = f->mol[i].bond[k].v2;
            t->coeff[1] = f->mol[i].bond[k].v1/2.0;
            t->coeff[2] = f->mol[i].bond[k].v3/3.0;
            t->coeff[3] = f->mol[i].bond[k].v4/4.0;
            break;
          default:
            msg_warn_f("bond potential \"%s\" not converted",
              dlpoly_prm_bond_pot_sym(f->mol[i].bond[k].pot,0));
          }
        /* potential ID */
        b->type = lammps_conv_pot_find(t,qt);
        queue_add(qb,b);
        }
      /* atom-ID offset */
      for (k=0; k<f->mol[i].n_atoms; k++)
        id += f->mol[i].atom[k].n_atoms;
      }
  /* convert potentials */
  lammps_conv_pot_save(&(d->bond),&(d->n_bonds),qb);
  lammps_conv_pot_save(&(d->pot_bond),&(d->n_bond_types),qt);
  /* clean memory */
  queue_free(qt);
  queue_free(qb);
  }

/* convert DLPOLY data to LAMMPS format
 
   c - DL_POLY structure data
   f - DL_POLY force-field data */
struct lammps_dat *lammps_from_dlpoly(struct dlpoly_cfg *c,
  struct dlpoly_fld *f) {
  double side[3],cosa[3],lxyz[3],tilt[3];
  unsigned i,j,k,l,m,n,t_id,m_id,a_id;
  short found;
  struct lammps_dat *d;
  struct lammps_type *t;
  struct queue *q;
  struct qdata *p;
  q = queue_alloc();
  d = lammps_dat_new();
  /* header */
  sprintf(d->header,"%s",c->header);
  /* simulation box */
  if (c->pbc && c->cell) {
    /* side lengths and cell angles */
    for (i=0; i<3; i++) {
      side[i] = vec_fnorm(c->cell[i],3);
      cosa[i] = vec_fprod_scalar(c->cell[(i+1)%3],c->cell[(i+2)%3],3);
      }
    /* size in x,y,z directions and tilting factors */
    lxyz[0] = side[0];
    tilt[0] = side[1]*cosa[2];
    tilt[1] = side[2]*cosa[1];
    lxyz[1] = sqrt((side[1]*side[1])-(tilt[0]*tilt[0]));
    tilt[2] = ((side[1]*side[2]*cosa[0])-(tilt[0]*tilt[1]))/lxyz[1];
    lxyz[2] = sqrt((side[2]*side[2])-(tilt[1]*tilt[1])-(tilt[2]*tilt[2]));
    /* save to LAMMPS format */
    d->box = lammps_box_new();
    for (i=0; i<3; i++) {
      d->box->pbc[i] = 1;
      d->box->min[i] = -0.5*lxyz[i];
      d->box->max[i] = 0.5*lxyz[i];
      d->box->tilt[i] = tilt[i];
      for (j=0; j<3; j++)
        d->box->vector[i][j] = c->cell[i][j];
      }
    }
  /* atoms */
  d->n_atoms = c->n_atoms;
  d->atom = lammps_atom_new(c->n_atoms);
  d->atom_style = LAMMPS_ATOM_FULL;
  a_id = 0;
  m_id = 0;
  for (i=0; i<f->n_molecules; i++) {
    for (j=0; j<f->mol[i].n_molecules; j++) {
      for (k=0; k<f->mol[i].n_atoms; k++) {
        for (l=0; l<f->mol[i].atom[k].n_atoms; l++) {
          /* type */
          n = atom_num_mass(f->mol[i].atom[k].mass);
          for (p=q->data, t_id=0, found=0; p; p=p->q_next, t_id++) {
            t = p->q_data;
            if (t->num==n && t->mass==f->mol[i].atom[k].mass &&
                str_compare(t->name,f->mol[i].atom[k].name)) {
              found = 1;
              break;
              }
            }
          if (!found) {
            t = lammps_type_new(1);
            sprintf(t->name,"%s",f->mol[i].atom[k].name);
            t->num = n;
            t->mass = f->mol[i].atom[k].mass;
            queue_add(q,t);
            }
          d->atom[a_id].type = t_id;
          /* molecular ID */
          d->atom[a_id].mol = m_id;
          /* charge */
          d->atom[a_id].charge = f->mol[i].atom[k].charge;
          /* coordinates */
          d->atom[a_id].crd = vec_falloc(3);
          for (m=0; m<3; m++)
            d->atom[a_id].crd[m] = c->crd[a_id][m];
          /* velocities */
          if (c->data && c->vel) {
            d->atom[a_id].vel = vec_falloc(3);
            for (m=0; m<3; m++)
              d->atom[a_id].vel[m] = c->vel[a_id][m];
            }
          a_id++;
          }
        }
      m_id++;
      }
    }
  /* atom types */
  d->n_atom_types = q->num;
  d->type = lammps_type_new(d->n_atom_types);
  for (i=0; i<d->n_atom_types; i++) {
    t = queue_get(q);
    sprintf(d->type[i].name,"%s",t->name);
    d->type[i].num = t->num;
    d->type[i].mass = t->mass;
    free(t);
    }
  /* bonding potentials */
  lammps_from_dlpoly_bonds(d,f);
  lammps_from_dlpoly_angles(d,f);
  lammps_from_dlpoly_dihedrals(d,f);
  lammps_from_dlpoly_pairs(d,f);
  /* clean memory */
  queue_free(q);
  return(d);
  }

/* -------------------------------------------------------------------------- */
