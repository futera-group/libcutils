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

#ifndef ZF_LIB_PRG_CP2K_H
#define ZF_LIB_PRG_CP2K_H

#include <stdio.h>
#include <mol/cell.h>
#include <qmc/basis.h>

/* reading flags */
#define READ_ALL          511  /* read all data */
#define READ_HEADER         1  /* header (basic system info) */
#define READ_ATOMS          2  /* atomic coordinates */
#define READ_KINDS          4  /* atomic kinds and basis set */
#define READ_OVERLAP        8  /* overlap matrix */
#define READ_ENERGY        16  /* energy */
#define READ_CHARGES       32  /* atomic charges */
#define READ_STATES        64  /* electronic states */
#define READ_ECOUPL       128  /* electronic coupling */
#define READ_ECOUPL_INFO  256  /* electronic coupling setup */

/* run types */
#define RUN_TYPE_BAND   1  /* band methods */
#define RUN_TYPE_BSSE   2  /* basis set superposition error */
#define RUN_TYPE_COPT   3  /* cell optimization */
#define RUN_TYPE_DEBG   4  /* debug analysis */
#define RUN_TYPE_DRVR   5  /* i-pi driver mode */
#define RUN_TYPE_EFST   6  /* ehrenfest dynamics */
#define RUN_TYPE_SPCT   7  /* electronic spectra */
#define RUN_TYPE_WOPT   8  /* energy calculation */
#define RUN_TYPE_EFRC   9  /* energy and forces */
#define RUN_TYPE_GOPT  10  /* geometry optimization */
#define RUN_TYPE_LINR  11  /* linear response */
#define RUN_TYPE_MNCR  12  /* monte carlo */
#define RUN_TYPE_MDYN  13  /* molecular dynamics */
#define RUN_TYPE_NONE  14  /* no task */
#define RUN_TYPE_NMOD  15  /* vibrational analysis */
#define RUN_TYPE_PINT  16  /* path integral */
#define RUN_TYPE_RTPR  17  /* real time propagation */
#define RUN_TYPE_TAMC  18  /* temperature accelerated monte carlo */
#define RUN_TYPE_TRMC  19  /* tree monte carlo */

/* output level */
#define OUT_TYPE_DEBG   1  /* debug - everything is written out */
#define OUT_TYPE_HIGH   2  /* high - lots of output */
#define OUT_TYPE_LOWX   3  /* low - litte output */
#define OUT_TYPE_MEDM   4  /* medium - quite some output */
#define OUT_TYPE_SLNT   5  /* silent - almost no output */

/* atomic charges */
#define CHRG_EFF        0  /* effective */
#define CHRG_MULL       1  /* Mulliken */
#define CHRG_LOWD       2  /* Lowdin */
#define CHRG_HIRSH      3  /* Hirshfeld */

/* cp2k atom data struct */
struct cp2k_atom {
  unsigned num;                    /* atomic number */
  unsigned kind;                   /* atom kind ID */
  double charge[4];                /* effective charge */
  double mass;                     /* atomic mass */
  double crd[3];                   /* coordinates [A] */
  };

/* cp2k atom kind data struct */
struct cp2k_kind {
  unsigned n_atoms;                /* number of atoms */
  unsigned n_shell_sets;           /* number of shell sets */
  unsigned n_shells;               /* number of shells */
  unsigned n_prim_fce;             /* number of primitive cartesian functions */
  unsigned n_cart_fce;             /* number of cartesian basis functions */
  unsigned n_sphr_fce;             /* number of spherical basis functions */
  short norm_type;                 /* orbital norm type */
  struct cp2k_shell *shell;        /* array of orbital shells */
  };

/* cp2k orbital shell data struct */
struct cp2k_shell {
  unsigned set;                    /* shell set ID */
  unsigned type;                   /* angular quantum number */
  unsigned n_cart_fce;             /* number of cartesian AOs */
  unsigned n_sphr_fce;             /* number of corresponding spher. fce */
  struct cp2k_ao *ao;              /* array of atomic orbitals */
  };

/* cp2k atom orbital data struct */
struct cp2k_ao {
  int qnum[3];                     /* principal/angular/magnetic quant. no. */
  unsigned n_prim_fce;             /* number of primitive functions */
  double *exp;                     /* primitive function exponents */
  double *coeff;                   /* primitive function coefficients */
  };

/* cp2k state data struct */
struct cp2k_state {
  double energy;                   /* state energy (KS eigenvalue) */
  double occup;                    /* state occupation */
  double *ao_pj;                   /* projection to AO basis set */
  };

/* cp2k electronic-coupling data block structure */
struct cp2k_ecpl_block {
  unsigned n_atoms;                /* number of atoms in the block */
  unsigned n_electrons;            /* number of electrons in the block */
  unsigned *atom_id;               /* array of atom IDs */
  unsigned *n_states;              /* number of states */
  double **energy;                 /* array of state energies */
  double ***ao_pj;                 /* projection to AO basis set */
  double ****coupling;             /* couplings (spin,block2,state1,state2) */
  };

/* cp2k electronic-coupling data structure */
struct cp2k_ecpl {
  unsigned n_blocks;               /* number of atomic blocks */
  unsigned *n_states;              /* number of states */
  double **energy;                 /* array of state energies */
  struct cp2k_ecpl_block * block;  /* array of the blocks */
  };

/* cp2k k-point data structure */
struct cp2k_kpoint {
  double crd[3];                   /* k-space coordinates */
  double weight;                   /* integration weight */
  struct cp2k_state **state;       /* array of electonic states [alpha/beta] */
  };

/* cp2k data struct */
struct cp2k_dat {
  char *version;                   /* program version */
  char *project_name;              /* project name */
  short run_type;                  /* run type */
  short print_level;               /* print level */
  unsigned charge;                 /* net charge */
  unsigned multiplicity;           /* spin multiplicity */
  unsigned n_atom_kinds;           /* number of atom kinds */
  unsigned n_atoms;                /* number of atoms */
  unsigned n_shell_sets;           /* number of shell sets */
  unsigned n_shells;               /* number of shells */
  unsigned n_prim_fce;             /* number of primitive cartesian functions */
  unsigned n_cart_fce;             /* number of cartesian basis functions */
  unsigned n_sphr_fce;             /* number of spherical basis functions */
  unsigned n_bfce;                 /* number of basis functions (AOs) */
  unsigned n_ibfce;                /* number of independent basis functions */
  unsigned n_spins;                /* number of spin components */
  unsigned n_kpoints;              /* number of k-points */
  unsigned *n_electrons;           /* number of electrons [alpha/beta] */
  unsigned *n_states;              /* number of states (MOs) [alpha/beta] */
  double fermi;                    /* Fermi energy [alpha/beta] */
  double **ovrl;                   /* AO overlap matrix */
  double **tmat;                   /* transformation matrix */
  struct cell *cell;               /* cell data */
  struct cp2k_atom *atom;          /* array of atoms */
  struct cp2k_kind *kind;          /* array of atom kinds */
  struct cp2k_ecpl *ecpl;          /* electronic coupling data */
  struct cp2k_kpoint *kpoint;      /* array of k-points */
  };

/* allocate and initialize array of atoms */
struct cp2k_atom* cp2k_atom_new(unsigned);
/* free memory allocated for cp2k atom array */
struct cp2k_atom* cp2k_atom_free(struct cp2k_atom*);

/* allocate and initialize array of atom kinds */
struct cp2k_kind* cp2k_kind_new(unsigned);
/* free memory allocated for cp2k atom kind array */
struct cp2k_kind* cp2k_kind_free(struct cp2k_kind*);

/* allocate and initialize array of orbital shells */
struct cp2k_shell* cp2k_shell_new(unsigned);
/* free memory allocated for cp2k orbital shell array */
struct cp2k_shell* cp2k_shell_free(struct cp2k_shell*, unsigned);

/* allocate and initialize array of atomic orbitals */
struct cp2k_ao* cp2k_ao_new(unsigned);
/* free memory allocated for cp2k atomic orbital array */
struct cp2k_ao* cp2k_ao_free(struct cp2k_ao*, unsigned);

/* allocate and initialize array of state data structs */
struct cp2k_state *cp2k_state_new(unsigned);
/* free memory allocated for array of state data structs */
struct cp2k_state** cp2k_state_free(struct cp2k_state**, unsigned*, unsigned);

/* allocate memory for electronic-coupling data structure */
struct cp2k_ecpl* cp2k_ecpl_new(void);
/* free memory allocated for cp2k electronic-coupling data structure */
struct cp2k_ecpl* cp2k_ecpl_free(struct cp2k_ecpl*, unsigned);
/* allocate memory for array electronic-coupling data blocks */
struct cp2k_ecpl_block* cp2k_ecpl_block_new(unsigned);
/* free memory allocated for array of cp2k electronic-coupling data blocks */
struct cp2k_ecpl_block* cp2k_ecpl_block_free(struct cp2k_ecpl_block*, 
  unsigned, unsigned);

/* allocate and initialize array of k-point data structs */
struct cp2k_kpoint *cp2k_kpoint_new(unsigned);
/* free memory allocated for array of k-point data structs */
struct cp2k_kpoint* cp2k_kpoint_free(struct cp2k_kpoint*, 
  unsigned, unsigned*, unsigned);

/* allocate and initialize new cp2k data struct */
struct cp2k_dat* cp2k_dat_new(void);
/* free memory allocated for cp2k data struct */
struct cp2k_dat* cp2k_dat_free(struct cp2k_dat*);

/* return internal ID of run type */
short cp2k_run_type_id(char*);
/* convert internal run type ID to corresponding keyword */
void cp2k_run_type_name(char*, short);
/* return internal ID of output type */
short cp2k_print_level_id(char*);
/* convert internal output type ID to corresponding keyword */
void cp2k_print_level_name(char*, short);
/* identify quantum numbers of from the orbital type keyword */
void cp2k_orbital_type_id(char*, int*, int*, int*);

/* read data from cp2k output file */
void cp2k_log_read(struct cp2k_dat*, int, char*);
/* read header data from cp2k output file */
void cp2k_log_read_header(struct cp2k_dat*, FILE*);
/* read atomic data from cp2k output file */
void cp2k_log_read_atoms(struct cp2k_dat*, FILE*);
/* read atomic kinds from cp2k output file */
void cp2k_log_read_kinds(struct cp2k_dat*, FILE*);
/* read overlap matrix from cp2k output file */
void cp2k_log_read_overlap(struct cp2k_dat*, FILE*);
/* read energies from cp2k output file */
void cp2k_log_read_energy(struct cp2k_dat*, FILE*);
/* read calculated atomic charges from cp2k output file */
void cp2k_log_read_charges(struct cp2k_dat*, FILE*);
/* read electronic states from cp2k output file */
void cp2k_log_read_states(struct cp2k_dat*, FILE*);
/* read electronic states from cp2k binary data file */
void cp2k_log_read_states_bin(struct cp2k_dat*, char*);
/* read electronic couplings from cp2k output file */
void cp2k_log_read_ecoupl(struct cp2k_dat*, int, FILE*);
/* read electronic couplings from cp2k binary data file */
void cp2k_log_read_ecoupl_bin(struct cp2k_dat*, char*, short);
/* read localized POD wavefunctions from cp2k binary data file */
void cp2k_log_read_ecoupl_wfn(struct cp2k_dat*, char*, short);
/* read POD transformation matrix from cp2k binary data file */
void cp2k_log_read_ecoupl_tmx(struct cp2k_dat*, char*, short);

/* convert cp2k structure to generic molecular format */
struct gen_mol *cp2k_dat_conv_mol(struct cp2k_dat*);
/* convert cp2k basis set to generic format */
struct basis *cp2k_dat_conv_basis(struct cp2k_dat*);

/* write structure coordinates into file */
void cp2k_dat_write_geom(struct cp2k_dat*, char*);
/* write structure coordinates into open file */
void cp2k_dat_write_geom_f(struct cp2k_dat*, FILE*);
/* write basis set into file */
void cp2k_dat_write_basis(struct cp2k_dat*, char*);
/* write basis set into open file */
void cp2k_dat_write_basis_f(struct cp2k_dat*, FILE*);
/* write overlap matrix into file */
void cp2k_dat_write_overlap(struct cp2k_dat*, char*);
/* write overlap matrix into open file */
void cp2k_dat_write_overlap_f(struct cp2k_dat*, FILE*);
/* write electronic states into file */
void cp2k_dat_write_states(struct cp2k_dat*, char*);
/* write electronic states into open file */
void cp2k_dat_write_states_f(struct cp2k_dat*, FILE*); 

#endif
