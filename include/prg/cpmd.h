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

#ifndef ZF_LIB_PRG_CPMD_H
#define ZF_LIB_PRG_CPMD_H

#include <stdio.h>
#include <cmn/queue.h>
#include <mol/cell.h>

/* basis functions */
#define CPMD_BASIS_S       1       /* S orbital */
#define CPMD_BASIS_P       3       /* P orbital - x,y,z */
#define CPMD_BASIS_D       5       /* D orbital */
#define CPMD_BASIS_F       7       /* F orbital */

/* structure format */
#define CPMD_FRMT_XYZ      0       /* XYZ file format */  
#define CPMD_FRMT_CIF      1       /* CIF file format */

/* charge type */
#define CPMD_CHARGE_ESP    0       /* fitted electrostatic potential */
#define CPMD_CHARGE_MUL    1       /* Mulliken charge */
#define CPMD_CHARGE_LOW    2       /* Lowdin charge */

/* cpmd atomic data struct */
struct cpmd_atom {
  unsigned t_id;            /* atomic type id */
  double charge[3];         /* atomic charge */
  double coord[3];          /* atomic coordinates */
  double valence;           /* atomic valence */
  };

/* cpmd atom type */
struct cpmd_type {
  unsigned num;             /* atomic number */
  unsigned n_ao_ids;        /* number of orbital types */
  unsigned n_orbs;          /* total number of atomic orbitals */
  double charge;            /* valence charge without PP core */
  short *ao_id;             /* type of atom orbitals */
  };

/* cpmd state data struct */
struct cpmd_state {
  double energy;            /* state energy (KS eigenvalue) */
  double compln;            /* completness of projection */
  double occup;             /* state occupation */
  double *ao_pj;            /* projection to AO basis set */
  };

/* cpmd k point data struct */
struct cpmd_kpoint {
  double coord[3];          /* K point coordinate */
  double weight;            /* weight */
  };

/* cpmd system structure data */
struct cpmd_geom {
  double energy;            /* structure energy */
  double **coord;           /* structure coordinates */
  };

/* cpmd data struct */
struct cpmd_dat {
  unsigned n_atoms;           /* number of atoms */
  unsigned n_types;           /* number of types */
  unsigned n_states;          /* number of states */
  unsigned n_alpha_states;    /* number of alpha states */
  unsigned n_beta_states;     /* number of beta states */
  unsigned n_electrons;       /* number of electrons */
  unsigned n_orbitals;        /* number of atomic orbitals */
  unsigned n_geoms;           /* number of structures from geom. optimization */
  unsigned n_kpoints;         /* number of K points */
  unsigned kp_mesh[3];        /* k-point mesh */
  unsigned multiplicity;      /* spin multiplicity */
  short lsd_calc;             /* LSD calculation */
  double charge;              /* total system charge */
  double **ovrl;              /* AO overlap matrix */
  struct cell *cell;          /* simulation cell */
  struct cpmd_type *type;     /* type array */
  struct cpmd_state *state_a; /* alpha state array */
  struct cpmd_state *state_b; /* beta state array */
  struct cpmd_kpoint *kpoint; /* kpoint array */
  struct cpmd_atom *atom;     /* atom array */
  struct cpmd_geom *g_final;  /* final structure */
  struct cpmd_geom *g_opt;    /* structures from geometry optimization */
  };

/* allocate and initialize new cpmd data struct */
struct cpmd_dat *cpmd_dat_new(void);
/* free memory allocated for cpmd data struct */
void cpmd_dat_free(struct cpmd_dat*);

/* allocate array of atom data structs */
struct cpmd_atom *cpmd_atom_new(unsigned);
/* free memory allocated for array of atom data structs */
struct cpmd_atom *cpmd_atom_free(struct cpmd_atom*, unsigned);

/* allocate array of geom data structs */
struct cpmd_geom *cpmd_geom_new(unsigned, unsigned);
/* free memory allocated for array of geom data structs */
struct cpmd_geom *cpmd_geom_free(struct cpmd_geom*, unsigned, unsigned);

/* allocate array of k-point data structs */
struct cpmd_kpoint *cpmd_kpoint_new(unsigned);
/* free memory allocated for array of k-point data structs */
struct cpmd_kpoint *cpmd_kpoint_free(struct cpmd_kpoint*, unsigned);

/* allocate array of state data structs */
struct cpmd_state *cpmd_state_new(unsigned);
/* free memory allocated for array of state data structs */
struct cpmd_state *cpmd_state_free(struct cpmd_state*, unsigned);

/* allocate array of type data structs */
struct cpmd_type *cpmd_type_new(unsigned);
/* free memory allocated for array of type data structs */
struct cpmd_type *cpmd_type_free(struct cpmd_type*, unsigned);

/* return internal code of orbital type (s,p,d,f) */
short cpmd_bs_id(char*);
/* convert inter code of orbital type to its symbol (s,p,d,f) */
char* cpmd_bs_sym(short, unsigned);
/* create array with atom orbital types for specified type of atom */
void cpmd_bs_set(struct cpmd_dat*, struct queue*, unsigned);

/* read data from cpmd output file */
void cpmd_log_read(struct cpmd_dat*, char*);
/* read atomic basis set from cpmd log file */
void cpmd_log_read_basis(struct cpmd_dat*, char*);
/* read atomic basis for from cpmd log file */
void cpmd_log_read_basis_f(struct cpmd_dat*, FILE*);
/* read system structure from cpmd log file */
void cpmd_log_read_geom(struct cpmd_dat*, char*);
/* read system structure from open cpmd log file */
void cpmd_log_read_geom_f(struct cpmd_dat*, FILE*);
/* read charge population from cpmd log file */
void cpmd_log_read_pop(struct cpmd_dat*, char*);
/* read charge population from open cpmd log file */
void cpmd_log_read_pop_f(struct cpmd_dat*, FILE*);
/* read projection of states from cpmd log file */
void cpmd_log_read_proj(struct cpmd_dat*, char*);
/* read projection of states from open cpmd log file */
void cpmd_log_read_proj_f(struct cpmd_dat*, FILE*);
/* read system specification from cpmd log file */
void cpmd_log_read_spec(struct cpmd_dat*, char*);
/* read system specification from open cpmd log file */
void cpmd_log_read_spec_f(struct cpmd_dat*, FILE*);
/* read state data from cpmd log file */
void cpmd_log_read_states(struct cpmd_dat*, char*);
/* read state data from open cpmd log file */
void cpmd_log_read_states_f(struct cpmd_dat*, FILE*);

/* read overlap matrix of AO orbitals from cpmd overlap file */
void cpmd_ovrl_read_bin(struct cpmd_dat*, char*);
/* read overlap matrix of AO orbitals from text-formatted cpmd overlap file */
void cpmd_ovrl_read_txt(struct cpmd_dat*, char*);

/* read projection of states from cpmd wavefunction file */
void cpmd_wfn_read_bin(struct cpmd_dat*, char*);
/* read projection of states from cpmd text-formatted wavefunction file */
void cpmd_wfn_read_txt(struct cpmd_dat*, char*);

/* set number of types and type ID for each atom */
void cpmd_dat_set_types(struct cpmd_dat*, unsigned*);

/* write energy bands into open file */
void cpmd_dat_write_band(struct cpmd_dat*, char*);
/* write energy bands into open file */
void cpmd_dat_write_band_f(struct cpmd_dat*, FILE*);
/* write structure coordinates into file */
void cpmd_dat_write_geom(struct cpmd_dat*, long int, short, char*);
/* write structure coordinates into open file */
void cpmd_dat_write_geom_f(struct cpmd_dat*, long int, short, FILE*);
/* write state projection into file */
void cpmd_dat_write_proj(struct cpmd_dat*, char*);
/* write state projection into open file */
void cpmd_dat_write_proj_f(struct cpmd_dat*, FILE*);

/* convert cpmd structure to generic molecular format */
struct gen_mol *cpmd_dat_conv_mol(struct cpmd_dat*, long);

#endif
