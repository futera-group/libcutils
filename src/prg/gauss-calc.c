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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/vector.h>
#include <qmc/basis.h>
#include <qmc/gto.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* calculate full overlap matrix over molecular orbitals
 
   g - the gaussian data struct
   t - spin
   s - the matrix */
void gauss_calc_mat_s(struct gauss_dat *g, short t, double **s) {
  unsigned i,j,k,l;
  double **r;
  struct gauss_mo *m;
  m = (t==GAUSS_SPIN_A ? g->mo_a : g->mo_b);
  /* AO overlap matrix */
  r = mat_falloc(g->bs->n_bfce,g->bs->n_bfce);
  gto_int_1e_s_mat(g->bs,r);
  /* MO overlap matrix */
  for (i=0; i<g->bs->n_ibfce; i++) {
    for (j=0; j<g->bs->n_ibfce; j++) {
      s[i][j] = 0.0;
      for (k=0; k<g->bs->n_bfce; k++)
        for (l=0; l<g->bs->n_bfce; l++)
          s[i][j] += (r[k][l]*m[i].coeff[k]*m[j].coeff[l]);
      }
    }
  /* clean memory */
  mat_ffree(r,g->bs->n_bfce);
  }

/* -------------------------------------------------------------------------- */

/* canonical orthogonalization of molecular orbitals
 
   g - the gaussian data struct
   t - spin */
void gauss_calc_mo_orth_cn(struct gauss_dat *g, short t) {
  unsigned i,j,k,n;
  double **s,**u,*e,**c;
  struct gauss_mo *m;
  if (g->bs->n_ibfce!=g->bs->n_bfce)
    msg_error("canonical orthogonalization of non-square MO matrix",1);
  n = g->bs->n_bfce;
  m = (t==GAUSS_SPIN_A ? g->mo_a : (g->mo_b ? g->mo_b : g->mo_a));
  /* overlap matrix */
  s = mat_falloc(n,n);
  gauss_calc_mat_s(g,t,s);
  /* diagonalization */
  u = mat_falloc(n,n);
  e = vec_falloc(n);
  mat_diag(s,u,e,n,0);
  /* orghogonalization */
  c = mat_falloc(n,n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c[i][j] = 0.0;
      for (k=0; k<n; k++)
        c[i][j] += u[k][i]*m[k].coeff[j];
      c[i][j] /= sqrt(e[i]);
      }
  for (i=0; i<g->bs->n_ibfce; i++)
    for (j=0; j<g->bs->n_bfce; j++)
      m[i].coeff[j] = c[i][j];
  /* clean memory */
  mat_ffree(c,n);
  mat_ffree(s,n);
  mat_ffree(u,n);
  vec_ffree(e);  
  }

/* symmetric Lowdin orthogonalization of molecular orbitals
 
   g - the gaussian data struct
   t - spin */
void gauss_calc_mo_orth_lw(struct gauss_dat *g, short t) {
  unsigned i,j,k,l,n;
  double **s,**u,*e,**c,v;
  struct gauss_mo *m;
  if (g->bs->n_ibfce!=g->bs->n_bfce)
    msg_error("canonical orthogonalization of non-square MO matrix",1);
  n = g->bs->n_bfce;
  m = (t==GAUSS_SPIN_A ? g->mo_a : (g->mo_b ? g->mo_b : g->mo_a));
  /* overlap matrix */
  s = mat_falloc(n,n);
  gauss_calc_mat_s(g,t,s);
  /* diagonalization */
  u = mat_falloc(n,n);
  e = vec_falloc(n);
  mat_diag(s,u,e,n,0);
  /* orghogonalization */
  c = mat_falloc(n,n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      c[i][j] = 0.0;
      for (k=0; k<n; k++) {
        v = 0.0;
        for (l=0; l<n; l++)
          v += (u[i][l]*u[k][l]/sqrt(e[l]));
        c[i][j] += (m[k].coeff[j]*v);
        }
      }
  for (i=0; i<g->bs->n_ibfce; i++)
    for (j=0; j<g->bs->n_bfce; j++)
      m[i].coeff[j] = c[i][j];
  /* clean memory */
  mat_ffree(c,n);
  mat_ffree(s,n);
  mat_ffree(u,n);
  vec_ffree(e);  
  }

/* orthogonalize molecular orbitals saved in the gausssian data struct
 
   g - the gaussian data struct
   t - orthogonalization algorithm */
void gauss_calc_mo_orth(struct gauss_dat *g, short t) {
  switch (t) {
    case GAUSS_MO_ORTH_CN:
      gauss_calc_mo_orth_cn(g,GAUSS_SPIN_A);
      gauss_calc_mo_orth_cn(g,GAUSS_SPIN_B);
      break;
    case GAUSS_MO_ORTH_LW:
      gauss_calc_mo_orth_lw(g,GAUSS_SPIN_A);
      gauss_calc_mo_orth_lw(g,GAUSS_SPIN_B);
      break;
    default:
      msg_error("invalid algoritm ID in gauss_calc_mo_orth",1);
      break;
    }
  }

/* sort molecular orbitals according to their energies
 
   g - the gaussian data struct */
void gauss_calc_mo_sort(struct gauss_dat *g) {
  double en,*cf;
  unsigned i,j,k,n;
  short chg;
  struct gauss_mo *m;
  n = g->bs->n_ibfce;
  if (n<2)
    return;
  /* alpha and beta spin */
  for (k=0; k<2; k++) {
    m = (k==0 ? g->mo_a : g->mo_b);
    if (!m)
      break;
    /* buble sort */
    for (i=1; i<n; i++) {
      chg = 0;
      for (j=n-1; j>=i; j--) {
        if (m[j-1].energy>m[j].energy) {
          /* swqp two MOs */
          en = m[j-1].energy;
          cf = m[j-1].coeff;
          m[j-1].energy = m[j].energy;
          m[j-1].coeff = m[j].coeff;
          m[j].energy = en;
          m[j].coeff = cf;
          }
        }
      if (chg)
        break;
      }
    }
  }

/* -------------------------------------------------------------------------- */
