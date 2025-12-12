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

#ifndef ZF_LIB_QMC_BASIS_H
#define ZF_LIB_QMC_BASIS_H

#include <stdio.h>

#define BASIS_DEF_GEN     0  /* general basis set */
#define BASIS_DEF_STO3G   1  /* minimal STO-3G basis set */

#define BASIS_SHELL_NUM  16  /* total number of shell types */
#define BASIS_SHELL_X     0  /* unknown shell type */
#define BASIS_SHELL_S     1  /* s-type shell */
#define BASIS_SHELL_SP    2  /* sp-type shell */
#define BASIS_SHELL_P     3  /* p-type shell */
#define BASIS_SHELL_Dp    4  /* pure d-type shell */
#define BASIS_SHELL_Dc    5  /* cartesian d-type shell */
#define BASIS_SHELL_Fp    6  /* pure f-type shell */
#define BASIS_SHELL_Fc    7  /* cartesian f-type shell */
#define BASIS_SHELL_Gp    8  /* pure g-type shell */
#define BASIS_SHELL_Gc    9  /* cartesian g-type shell */
#define BASIS_SHELL_Hp   10  /* pure h-type shell */
#define BASIS_SHELL_Hc   11  /* cartesian h-type shell */
#define BASIS_SHELL_Ip   12  /* pure i-type shell */
#define BASIS_SHELL_Ic   13  /* cartesian i-type shell */
#define BASIS_SHELL_Jp   14  /* pure j-type shell */
#define BASIS_SHELL_Jc   15  /* cartesian j-type shell */

/* basis set shell struct */
struct basis_shell {
  short type;                 /* shell type */
  unsigned bf;                /* basis function ID */
  unsigned n_bfce;            /* number of basis functions */
  unsigned n_prim;            /* number of primitive functions */
  double *exp;                /* exponents of primitive functions */
  double *cf1;                /* coefficients of primitive functions */
  double *cf2;                /* coefficients of primitive functions (SP) */
  };

/* basis set center struct */
struct basis_center {
  unsigned bf;                /* basis function ID */
  unsigned type;              /* atomic number */
  unsigned n_bfce;            /* number of basis functions */
  unsigned n_shells;          /* number of shells */
  double coord[3];            /* cartesian coordinates */
  struct basis_shell *shell;  /* shells for each type of functions */
  };

/* basis set struct */
struct basis {
  unsigned n_bfce;             /* number of basis functions */
  unsigned n_ibfce;            /* number of independent functions */
  unsigned n_centers;          /* number of basis set centers */
  unsigned n_cont_shells;      /* number of contracted shells */
  unsigned n_prim_shells;      /* number of primitive shells */
  unsigned n_pure_d;           /* number of pure d shells */
  unsigned n_pure_f;           /* number of pure f shells */
  unsigned max_ang_mom;        /* highest angular momentum */
  unsigned max_bfce_cont;      /* largest degree or contraction */
  struct basis_center *center; /* basis set centers */
  };

/* allocate memory for basis set struct */
struct basis* basis_new(void);
/* allocate memory for basis set shell struct */
struct basis_shell* basis_shell_new(unsigned);
/* allocate memory for basis set center struct */
struct basis_center* basis_center_new(unsigned);

/* free memory allocated for basis set */
struct basis* basis_free(struct basis*);
/* free memory allocated for basis set shell struct */
struct basis_shell* basis_shell_free(struct basis_shell*, unsigned);
/* free memory allocated for basis set center struct */
struct basis_center* basis_center_free(struct basis_center*, unsigned);

/* copy data from one basis set shell data struct to another */
void basis_shell_copy(struct basis_shell*, struct basis_shell*);
/* create copy of array of basis set shell data structure */
struct basis_shell* basis_shell_copy_new(struct basis_shell*, unsigned);
/* copy data from one basis set center data struct to another */
void basis_center_copy(struct basis_center*, struct basis_center*);
/* copy array of basis set center data structures */
struct basis_center* basis_center_copy_new(struct basis_center*, unsigned);

/* merge two basis sets to one basis set data struct */
struct basis* basis_merge(struct basis*, struct basis*);

/* return number of primitive basis set functions */
unsigned basis_bfce_nprim(struct basis*);
/* locate basis function according to its ID and return points to it */
void basis_bfce_pointer(struct basis*, unsigned, struct basis_center**,
  struct basis_shell**, unsigned*);

/* return ID of the shell which name is given */
short basis_shell_id(char*, short);
/* convert shell id to string name */
char* basis_shell_name(short);
/* convert shell id to string function name */
char* basis_shell_name_fce(short, unsigned);
/* return angular momentum of specified shell */
unsigned basis_shell_ang_mom(short);
/* return number of basis functions in the shell */
unsigned basis_shell_nbfce(short);
/* return number of shells in the basis set */
unsigned basis_shell_num(struct basis*);

/* define specified basis set */
void basis_def(struct basis*, struct basis*, short);
/* define user-defined basis set */
void basis_def_gen(struct basis*, struct basis*);
/* define minimal STO-3G basis set */
void basis_def_sto3g(struct basis*);
/* define global parameters of basis set */
void basis_def_parm(struct basis*);
/* convert basis set type ID to string name */
char* basis_def_name(short);
/* return ID of basis set type */
short basis_def_id(char*);

/* print basis set info */
void basis_print(struct basis*, short, char*);
/* print basis set info to file*/
void basis_fprint(struct basis*, short, FILE*);

/* write basis set in C format to file */
void basis_c_write(struct basis*, char*);

/* read basis set in CPMD format from file */
void basis_cpmd_read(struct basis*, char*);
/* write basis set in CPMD format to file */
void basis_cpmd_write(struct basis*, char*);

/* write basis set in gamess format to file */
void basis_gamess_write(struct basis*, char*);

/* read basis set in gaussian format from file */
void basis_gauss_read(struct basis*, char*);
/* read basis set in gaussian format from file */
void basis_gauss_read_f(struct basis*, FILE*);
/* write basis set in gaussian format to file */
void basis_gauss_write(struct basis*, char*);

/* write basis set in molpro format to file */
void basis_molpro_write(struct basis*, char*);

/* read basis set in NBO format from file */
void basis_nbo_read(struct basis*, char*);
/* write basis set in NBO format to file */
void basis_nbo_write(struct basis*, char*);

/* write basis set in turbomole format to file */
void basis_turbo_write(struct basis*, char*);

#endif
