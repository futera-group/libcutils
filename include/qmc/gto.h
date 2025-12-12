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

#ifndef ZF_LIB_QMC_GTO_H
#define ZF_LIB_QMC_GTO_H

#include "qmc/basis.h"

/* spherical -> cartesian GTO transformation */
#define GTO_LC_MAX      6    /* maximum number of terms */
#define GTO_MAX_AMOM    4    /* maximal allowed angular momentum */
#define GTO_MAX_I2EXP   81   /* 2e integral expansion maximum */

/* GTO pair PQ term */
struct gto_pair_pq {
  unsigned id[3];            /* term ID */
  double coeff;              /* expansion coefficient */
  };

/* GTO pair data struct for integrals */
struct gto_pair {
  unsigned aa1[3],aa2[3];    /* angular momenta (vector) */
  unsigned am1,am2;          /* angular momenta (modulus) */
  double e1,e2,e12;          /* exponents and their sum */
  double sigma;              /* sigma factor = 1/(2*e12) */
  double pref;               /* integral prefactor */
  double r12[3];             /* intershell vector */
  double crd[3];             /* pair center coordinates */
  unsigned n_pq;             /* number of PQ terms */
  struct gto_pair_pq pq[GTO_MAX_I2EXP];
  };

/* allocate memory for new GTO pair array */
struct gto_pair* gto_pair_new(unsigned);
/* free memory allocated for GTO pair array */
struct gto_pair* gto_pair_free(struct gto_pair*);

/* calculate and GTO pair integral data and select only important terms */
void gto_pair_eval_all(struct gto_pair*, struct basis*, double*, struct basis*,
  double*, long unsigned*, long unsigned, double);
/* evaluated one pair GTO data */
void gto_pair_eval_one(struct gto_pair*, struct basis_shell*, double*,
  unsigned, unsigned *, struct basis_shell*, double*, unsigned, unsigned*,
  double*, double, double, long unsigned, long unsigned*,
  long unsigned, double);
/* create copy of GTO pair with swapped data */
void gto_pair_swap(struct gto_pair*, struct gto_pair*);

/* get angular vector of specific cartesian gaussian type orbital */
void gto_ang_vec(short, short, unsigned*);
/* get name of specified angular momemntum */
char gto_ang_name(unsigned);

/* transform spherical GTO to linear combination of Cartesian GTO */
void gto_pure_to_cart(short, short, unsigned*, double*, unsigned*);

/* calculate normalization constant for Cartesian GTO */
double gto_norm(double, unsigned*);

/* evaluate auxiliary f(l,m,a,b) function */
double gto_afce_f(unsigned, unsigned, unsigned, double, double);
/* evaluate auxiliary A(l1,l2,A,B,C,g) function */
double gto_afce_a(unsigned, unsigned, unsigned, unsigned, unsigned,
  double, double, double, double);

/* calculate values of all basis functions at specific point */
double gto_bfce_val(struct basis_center*, struct basis_shell*,
  unsigned, double*);
/* calculate values of all basis functions at specific point */
void gto_bfce_val_all(struct basis*, double*, double*);

/* calculate overlap matrix */
void gto_int_1e_s_mat(struct basis*, double**);
/* calculate overlap matrix between two basis sets */
void gto_int_1e_s_mat_b12(struct basis*, struct basis*, double**);
/* calculate overlap integral between two gaussian basis functions */
double gto_int_1e_s_bf(struct basis_shell*, double*, unsigned,
  struct basis_shell*, double*, unsigned, double*, double);
/* calculate overlap integral between two core GTO functions */
double gto_int_1e_s_core(double, unsigned*, double*, double,
  unsigned*, double*);

/* calculate kinetic energy matrix */
void gto_int_1e_t_mat(struct basis*, double**);
/* calculate kinetic energy integral between two gaussian basis functions */
double gto_int_1e_t_bf(struct basis_shell*, double*, unsigned,
  struct basis_shell*, double*, unsigned, double*, double);
/* calculate kinetic energy integral between two core GTO functions */
double gto_int_1e_t_core(double, unsigned*, double*, double,
  unsigned*, double*);

/* calculate potential energy matrix */
void gto_int_1e_v_mat(struct basis*, double**, double*, unsigned, double**);
/* calculate potential energy integral between two gaussian basis functions */
double gto_int_1e_v_bf(struct basis_shell*, double*, unsigned,
  struct basis_shell*, double*, unsigned, double*, double,
  double*, double);
/* calculate potential energy integral between all gaussian basis functions */
void gto_int_1e_v_bf_all(struct basis*, double*, double**);
/* calculate potential energy integral between two core GTO functions */
double gto_int_1e_v_core(unsigned*, double *, unsigned*, double*,
  double, double*, double);

/* calculate HF 2e integral matrix */
void gto_int_2e_g_mat(struct basis*, double**, double**,
  double, double, double);
/* calculate one element of HF 2e integral matrix */
void gto_int_2e_g_elm(struct basis*, double**, struct basis_shell*,
  double*, unsigned, struct basis_shell*, double*, unsigned,
  double*, double*, double, double, double);
/* coulombic and exchange two electron integrals for one MO pair */
void gto_calc_2e_jk(struct basis*, double*, struct basis*, double*, 
  double*, double*, double, double, double);
/* evaluate one (ab|cd) 2e integral from 4 contracted GTO shells */
double gto_int_2e_e_bf(struct basis_shell*, double*, unsigned,
  struct basis_shell*, double*, unsigned, struct basis_shell*, double*,
  unsigned, struct basis_shell*, double*, unsigned, double*, double, 
  double*, double, double, double, double);
/* evaluate one (ab|cd) 2e integral from 4 contracted GTO shells */
double gto_int_2e_e_bf_r(
  struct basis_shell*, double*, unsigned, struct basis_shell*, double*,
  unsigned, struct basis_shell*, double*, unsigned, struct basis_shell*,
  double*, unsigned, double up, double uu, double zr);
/* calculate argument for Boys function */
void gto_int_2e_e_parm(struct gto_pair*, struct gto_pair*,
  double*, double*, double*);
/* evaluate one two electron integral from pre-calculated data */
double gto_int_2e_e_core(struct gto_pair*, struct gto_pair*,
  double*, double, double, double, double);
/* evaluate set of [r]^m auxiliary integrals from [0]^m integrals */
short gto_int_2e_rm(
  double rm[2*(GTO_MAX_AMOM+1)+1][2*(GTO_MAX_AMOM+1)+1]
           [2*(GTO_MAX_AMOM+1)+1][4*(GTO_MAX_AMOM+1)+1],
  unsigned, double, double, double, double*);
/* expand GTO pair bra or ket into set of p-bras of q-kets */
void gto_int_2e_pq(struct gto_pair*, double);
/* get name of the 2e integral class */
char *gto_int_2e_name(unsigned, unsigned, unsigned, unsigned);
/* auxiliary function for printing value of 2e integrals */
void gto_int_2e_print(struct basis_shell*, double*, unsigned, 
  struct basis_shell*, double*, unsigned, struct basis_shell*, double*, 
  unsigned, struct basis_shell*, double*, unsigned, double, double, 
  double, FILE*);
/* print out array with values of all symmetrically unique 2e AO integrals */
void gto_int_2e_print_all(struct basis*, double, double, double, FILE*);

#endif
