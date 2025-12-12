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
#include <cmn/stack.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* representation of GTO bra / ket pair */
struct gto_pq {
  unsigned a1[3],a2[3],a3[3]; /* angular vectors a,b,p */
  unsigned A1,A2,A3;          /* angular momentums a,b,p */
  unsigned c1,c2,c3;          /* expansion coefficients a',b',p' */
  double cf;                  /* coefficient */
  };

/* -------------------------------------------------------------------------- */

/* allocate memory for one PQ term data struct
 
   a1,a2,a3 - angular vectors a,b,p */
struct gto_pq* gto_int_2e_pq_new(unsigned *a, unsigned *b, unsigned *p) {
  unsigned i;
  struct gto_pq* t = NULL;
  /* allocate memory */
  t = (struct gto_pq*)malloc(sizeof(struct gto_pq));
  if (!t) 
    msg_error("cannot allocate memory for PQ term",1);
  /* save data */
  t->A1 = t->A2 = t->A3 = 0;
  t->c1 = t->c2 = t->c3 = 0;
  for (i=0; i<3; i++) {
    t->a1[i] = a[i];
    t->A1 += t->a1[i];
    t->a2[i] = b[i];
    t->A2 += t->a2[i];
    t->a3[i] =p[i];
    t->A3 += t->a3[i];
    }
  t->cf = 1.0;
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* expand GTO pair bra or ket according to the angular vector

   p - pointer to GTO pair data struct
   s - storage for the expansion
   t - term that is to be expanded
   c - expansion coordinate
   n - ID of angular vector */
void gto_int_2e_pq_a(struct gto_pair *p, struct stack *s, struct gto_pq *t,
  unsigned c, short n) {
  struct gto_pq *r;
  /* first term */
  if (t->a3[c]) {
    r = gto_int_2e_pq_new(t->a1,t->a2,t->a3);
    switch (n) {
      case 1:
        r->a1[c] -= 1;
        r->A1 -= 1;
        break;
      case 2:
        r->a2[c] -= 1;
        r->A2 -= 1;
        break;
      default:
        msg_error_f("invalid angular vector ID (%d) in GTO pair expansion",1,n);
        break;
      }
    r->a3[c] -= 1;
    r->A3 -= 1;
    r->c1 = t->c1;
    r->c2 = t->c2;
    r->c3 = t->c3;
    r->cf = t->cf*t->a3[c];
    stack_add(s,r);
    }
  /* second term */
  if (fabs(p->r12[c])>0.0) {
    r = gto_int_2e_pq_new(t->a1,t->a2,t->a3);
    switch (n) {
      case 1:
        r->a1[c] -= 1;
        r->A1 -= 1;
        r->c1 = t->c1;
        r->c2 = t->c2+1;
        r->cf = -t->cf;
        break;
      case 2:
        r->a2[c] -= 1;
        r->A2 -= 1;
        r->c1 = t->c1+1;
        r->c2 = t->c2;
        r->cf = t->cf;
        break;
      default:
        msg_error_f("invalid angular vector ID (%d) in GTO pair expansion",1,n);
        break; }
    r->c3 = t->c3+1;
    r->cf *= p->r12[c];
    stack_add(s,r);
  }
  /* third term */
  r = gto_int_2e_pq_new(t->a1,t->a2,t->a3);
  switch (n) {
    case 1:
      r->a1[c] -= 1;
      r->A1 -= 1;
      break;
    case 2:
      r->a2[c] -= 1;
      r->A2 -= 1;
      break;
    default:
      msg_error_f("invalid angular vector ID (%d) in GTO pair expansion",1,n);
      break;
    }
  r->a3[c] += 1;
  r->A3 += 1;
  r->c1 = t->c1;
  r->c2 = t->c2;
  r->c3 = t->c3+1;
  r->cf = t->cf;
  stack_add(s,r);
  }

/* expand GTO pair bra or ket into set of p-bras of q-kets
 
   p  - pointer to GTO pair data struct 
   zr - criterion for neglecting the contribution */
void gto_int_2e_pq(struct gto_pair *p, double zr) {
  unsigned i,a0[3] = {0,0,0};
  struct gto_pq *t;
  struct stack *s;
  s = stack_alloc();
  /* initial bra / ket */
  t = gto_int_2e_pq_new(p->aa1,p->aa2,a0);
  stack_add(s,t);
  /* expansion */
  p->n_pq = 0;
  while (s->num) {
    t = stack_get(s);
    /* expand according to the first vector */
    if (t->A1) {
      for (i=0; i<3; i++)
        if (t->a1[i]) {
          gto_int_2e_pq_a(p,s,t,i,1);
          break;
          }
      }
    /* expand according to the second vector */
    else if (t->A2) {
      for (i=0; i<3; i++)
        if (t->a2[i]) {
          gto_int_2e_pq_a(p,s,t,i,2);
          break;
          }
      }
    /* save final RM term */
    else {
      if (fabs(t->cf)>0.0) {
        /* angular vector */
        for (i=0; i<3; i++)
          p->pq[p->n_pq].id[i] = t->a3[i];
        /* coefficient */
        p->pq[p->n_pq].coeff = t->cf;
        p->pq[p->n_pq].coeff *= pow(2.0*p->e1,t->c1);
        p->pq[p->n_pq].coeff *= pow(2.0*p->e2,t->c2);
        p->pq[p->n_pq].coeff /= pow(2.0*p->e12,t->c3);
        /* add to array */
        if (t->cf>zr) {
          p->n_pq++;
          if (p->n_pq==GTO_MAX_I2EXP)
            msg_error_f("expansion maximum (%d) of ((%d,%d,%d),(%d,%d,%d)| ket"
              " reached",1,GTO_MAX_I2EXP,
              p->aa1[0],p->aa1[1],p->aa1[2],
              p->aa2[0],p->aa2[1],p->aa2[2]);
          }
        }
      }
    /* clean memory */
    free(t);
    }
  /* clean memory */
  stack_free(s);
  }

/* -------------------------------------------------------------------------- */
