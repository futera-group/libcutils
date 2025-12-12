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

#ifndef ZF_LIB_QMC_PSEUDO_H
#define ZF_LIB_QMC_PSEUDO_H

/* pseudopotential types */
#define PSEUDO_NC_NUM      1  /* norm-conserving numeric */
                           
/* DFT functionals */      
#define PSEUDO_XC_LYP      1  /* Lee-Yang_Par */
#define PSEUDO_XC_B88      2  /* Becke (1988) */
                           
/* electron shells */      
#define PSEUDO_L_S         0  /* S shell - 2 e */
#define PSEUDO_L_P         1  /* P shell - 6 e */
#define PSEUDO_L_D         2  /* D shell - 10 e */
#define PSEUDO_L_F         3  /* F shell - 14 e */
#define PSEUDO_L_G         4  /* G shell - 18 e */

/* pseudo-potential data structure */
struct pseudo {
  unsigned atom_num;                /* atomic number */
  unsigned val_num;                 /* number of valence electrons */
  unsigned xc_num;                  /* exchange-correlation - number */
  double xc_slater;                 /* exchange-correlation - slater */
  short xc_lda_c;                   /* exchange-correlation - LDA corr */
  short xc_gc_e;                    /* exchange-correlation - GC exchange */
  short xc_gc_c;                    /* exchange-correlation - GC corr */
  short pp_type;                    /* pseudopotential type */
  char *time_stamp;                 /* time stamp */
  unsigned state_core;              /* number of core states */
  unsigned state_val;               /* number of valence states */
  unsigned mesh_num;                /* number of mesh points */
  double en_pot;                    /* pseudoatom total energy */
  double en_tot;                    /* pseudoatom total energy */
  unsigned *el_conf_n;              /* electron configuration - n */
  unsigned *el_conf_l;              /* electron configuration - l */
  double *el_conf_occ;              /* electron configuration - occupation */
  unsigned mt_num;                  /* MT PP - bumbef if shells */
  unsigned *mt_pp_n;                /* MT PP - shell */
  short *mt_pp_l;                   /* MT PP - shell type */
  double *mt_pp_r;                  /* MT PP - core radius */
  double *mt_pp_e;                  /* MT PP - energy */
  double **dat_pot;                 /* potential data */
  double **dat_wfce;                /* wavefunction data */
  double **dat_den;                 /* atom density data */
  };

/* allocate memory for pseudopotential data */
struct pseudo *pseudo_new(void);
/* free memory allocated for pseudopotential data */
void pseudo_free(struct pseudo*);

/* return ID of electronic shell */
short pseudo_shell_id(char*);
/* return name of electronic shell */
char *pseudo_shell_name(short);

/* return ID of exchange-correlation functional */
short pseudo_func_id(char*);
/* return name of exchange-correlation functional */
char *pseudo_func_name(short);

/* read pseudopotential data from file */
void pseudo_read(struct pseudo*, char*);

/* write pseudopotential data into the file */
void pseudo_write(struct pseudo*, char*);

#endif
