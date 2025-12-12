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

#ifndef ZF_LIB_PRG_GAUSS_H
#define ZF_LIB_PRG_GAUSS_H

#include <stdio.h>

/* symbolic constants */
#define GAUSS_PRG_UNK       0  /* unknown program code */
#define GAUSS_PRG_98       98  /* Gaussian 98 */
#define GAUSS_PRG_03        3  /* Gaussian 03 */
#define GAUSS_PRG_09        9  /* Gaussian 09 */
#define GAUSS_PRG_A01     101  /* Revision A.01 */
#define GAUSS_PRG_B01     201  /* Revision B.01 */
#define GAUSS_PRG_C01     301  /* Revision C.01 */

#define GAUSS_MTHD_UNK      0  /* unknown method */
#define GAUSS_MTHD_AM1      1  /* AM1 */
#define GAUSS_MTHD_PM3      2  /* PM3 */
#define GAUSS_MTHD_PM6      3  /* PM6 */
#define GAUSS_MTHD_CNDO     4  /* CNDO */
#define GAUSS_MTHD_INDO     5  /* INDO */
#define GAUSS_MTHD_MINDO    6  /* MINDO */
#define GAUSS_MTHD_ZINDO    7  /* ZINDO */
#define GAUSS_MTHD_HF      10  /* HF */
#define GAUSS_MTHD_CIS     20  /* CIS */
#define GAUSS_MTHD_CISd    21  /* CIS(D) */
#define GAUSS_MTHD_CISD    22  /* CISD */
#define GAUSS_MTHD_CID     23  /* CID */
#define GAUSS_MTHD_CCSD    24  /* CCSD */
#define GAUSS_MTHD_CCSDt   25  /* CCSD(T) */
#define GAUSS_MTHD_EOMCCSD 26  /* EOM-CCSD */
#define GAUSS_MTHD_MP2     27  /* MP2 */
#define GAUSS_MTHD_MP3     28  /* MP3 */
#define GAUSS_MTHD_DFT    100  /* DFT */
#define GAUSS_MTHD_TDDFT  101  /* TDDFT */

#define GAUSS_SPIN_A        0  /* alpha spin */
#define GAUSS_SPIN_B        1  /* beta spin */

#define GAUSS_CRD_ANY       0  /* input or standard orientation coordinates */
#define GAUSS_CRD_INP       1  /* input orientation coordinates */
#define GAUSS_CRD_STD       2  /* standard orientation coordinates */
#define GAUSS_CRD_ZMT       3  /* z-matrix orientation coordinates */
#define GAUSS_CRD_TRJ       4  /* trajectory type of coordinate record */

#define GAUSS_CH_NC         0  /* nuclear charge */
#define GAUSS_CH_NEFF       1  /* effective nuclear charge */
#define GAUSS_CH_MM         2  /* MM charge */
#define GAUSS_CH_MM_MOD     3  /* modified MM charge */
#define GAUSS_CH_MULL       4  /* Mulliken charge */

#define GAUSS_JOB_UNK       0  /* unknown type of job */
#define GAUSS_JOB_SP        1  /* single point */
#define GAUSS_JOB_FOPT      2  /* full optimization */
#define GAUSS_JOB_POPT      3  /* partial optimization */
#define GAUSS_JOB_FTS       4  /* full TS optimization */
#define GAUSS_JOB_PTS       5  /* partial TS optimization */
#define GAUSS_JOB_FSDD      6  /* full 2nd order saddle point opt */
#define GAUSS_JOB_PSDD      7  /* partial 2nd order saddle point opt */
#define GAUSS_JOB_FORCE     8  /* gradient calculation */
#define GAUSS_JOB_FREQ      9  /* vibrational frequencies */
#define GAUSS_JOB_SCAN     10  /* potential surface scan */
#define GAUSS_JOB_GUESS    11  /* MO generation only */
#define GAUSS_JOB_LST      12  /* linear synchronous transit */
#define GAUSS_JOB_STAB     13  /* test of SCF/KS stability */
#define GAUSS_JOB_RARCH    14  /* archive info generation */
#define GAUSS_JOB_MIXED    15  /* mixed method models */

#define GAUSS_BAS_BFCE      0  /* all basis functions */
#define GAUSS_BAS_SHELL     1  /* contracted shells */

#define GAUSS_MO_ORTH_CN    1  /* Canonical MO orthogonalization */
#define GAUSS_MO_ORTH_LW    2  /* Lowdin symmetrical MO orthogonalization */

#define GAUSS_IOP_MAX      48  /* total number of defined IOps */
#define GAUSS_IOP_UNK       0  /* unknown IOp */
#define GAUSS_IOP_01_038    1  /* entry control */
#define GAUSS_IOP_01_172    2  /* print internal coordinates */
#define GAUSS_IOP_02_012    3  /* crowding abort control */
#define GAUSS_IOP_02_017    4  /* symmetry distance criteria */
#define GAUSS_IOP_02_018    5  /* symmetry non-distance criteria */
#define GAUSS_IOP_02_040    6  /* initial structure saving */
#define GAUSS_IOP_03_005    7  /* basis set type */
#define GAUSS_IOP_03_006    8  /* number of gaussian basis functions */
#define GAUSS_IOP_03_007    9  /* diffuse and polarization functions */
#define GAUSS_IOP_03_011   10  /* 2e integral storage format */
#define GAUSS_IOP_03_016   11  /* pseudopotential */
#define GAUSS_IOP_03_017   12  /* specification of pseudopotential */
#define GAUSS_IOP_03_024   13  /* gaussian function printing (basis set) */
#define GAUSS_IOP_03_025   14  /* number of last 2e integral */
#define GAUSS_IOP_03_030   15  /* 2e integral symmetry */
#define GAUSS_IOP_03_041   16  /* semi-empirical methods */
#define GAUSS_IOP_03_074   17  /* exchange-correlation functional */
#define GAUSS_IOP_03_082   18  /* fitting density basis set for couloumb int */
#define GAUSS_IOP_03_116   19  /* type of SCF */
#define GAUSS_IOP_04_005   20  /* type of initial guess */
#define GAUSS_IOP_04_035   21  /* overlap matrix */
#define GAUSS_IOP_05_005   22  /* direct SCF control */
#define GAUSS_IOP_05_007   23  /* maximum number of iterations */
#define GAUSS_IOP_05_008   24  /* selection of the procedure of direct min */
#define GAUSS_IOP_05_013   25  /* action on convergence failure */
#define GAUSS_IOP_05_032   26  /* sleazy SCF */
#define GAUSS_IOP_05_035   27  /* orthonormal basis set */
#define GAUSS_IOP_05_038   28  /* varying integral cutoff during direct SCF */
#define GAUSS_IOP_06_007   29  /* MO printing */
#define GAUSS_IOP_06_008   30  /* density matrix printing */
#define GAUSS_IOP_06_009   31  /* full population matrix printing */
#define GAUSS_IOP_06_010   32  /* Gross orbital changes printing */
#define GAUSS_IOP_06_028   33  /* SCF density as current density */
#define GAUSS_IOP_08_006   34  /* bucket selection */
#define GAUSS_IOP_08_009   35  /* debug control */
#define GAUSS_IOP_08_010   36  /* window selection */
#define GAUSS_IOP_08_068   37  /* EOM-CCSD setting */
#define GAUSS_IOP_08_108   38  /* EOM-CCSD setting */
#define GAUSS_IOP_09_005   39  /* post HF method */
#define GAUSS_IOP_09_023   40  /* orbital localization */
#define GAUSS_IOP_09_041   41  /* number of states for Davidson diag. */
#define GAUSS_IOP_09_042   42  /* method and matrix block for diag. */
#define GAUSS_IOP_09_044   43  /* density matrix filling */
#define GAUSS_IOP_09_049   44  /* initial guess vectors in L914 */
#define GAUSS_IOP_09_068   45  /* initial guess for EOM-CC */
#define GAUSS_IOP_99_005   46  /* checkpoint file handling */
#define GAUSS_IOP_99_009   47  /* dipole and el. field archiving */

/* gaussian atomic data */
struct gauss_atom {
  double coord[3];                  /* cartesian coordinates */
  double coord_fix[3];              /* constrained coordinates */
  double grad[3];                   /* cartesian gradient */
  unsigned num;                     /* atomic numbers */
  double charge[5];                 /* atomic charge */
  double mass;                      /* real atomic weights */
  char type_name[2][25];            /* atom types - standard and modified */
  unsigned type_id[2];              /* atom type IDs */
  unsigned weight;                  /* integer atomic weights */
  unsigned res_info;                /* atom residue info */
  unsigned res_num;                 /* atom residue numbers */
  unsigned frag_info;               /* atom fragment info */
  unsigned nuc_spin;                /* nuclear spins */
  double nuc_qmom;                  /* nuclear Q momentum */
  double nuc_gfac;                  /* nuclear G factor */
  };

/* gaussian molecular orbital data */
struct gauss_mo {
  double energy;                    /* orbital energy */
  double *coeff;                    /* expansion to basis functions */
  };

/* gaussian data struct - states */
struct gauss_state {
  double energy;             /* energy of the state [eV] */
  double *dipole;            /* ground-to-excited-state transition dipole [au] */
  double **dipole_ee;        /* excited-to-excited-state transition dipole [au] */
  double lambda;             /* corresponding wavelength [nm] */
  double strength;           /* oscillatory strength */
  double spin_s2;            /* spin <s**2> */
  double rot_vel;            /* rotatory strength (velocity) */
  double rot_len;            /* rotatory strength (length) */
  };

/* gaussian data struct - structures */
struct gauss_geom {
  double energy;             /* structure energies */
  double **coord;            /* structure coordinates */
  };

/* gaussian data struct - nmr shifts */
struct gauss_nmr {
  double isotropic;          /* isotropic moment */
  double anisotropy;         /* anisotropic moment */
  double **tensor;           /* shield tensor */
  double *eigen;             /* tensor eigenvalues */
  };

/* gaussian data struct - vibrational frequencies */
struct gauss_freq {
  double freq;               /* frequency [cm-1] */
  double **displ;            /* coordinate displacement */
  double *mu_e;              /* electric dipole transition moment */
  double *mu_m;              /* magnetic dipole transition moment */
  double *mu_d;              /* dipole derivative */
  double *pol;               /* contribution to vibrational polarizability */
  double rmass;              /* reduced mass [AMU] */
  double fconst;             /* force constant */
  double ir;                 /* IR intensity [KM/mol] */
  double raman;              /* Raman scattering activity [A^4/AMU] */
  double depol_p;            /* depolarization ratio (plane inc. light) */
  double depol_u;            /* depolarization ratio (unpolar. inc. light) */
  double strength_dip;       /* dipole strength [10^(-40) esu^2 cm^2] */
  double strength_rot;       /* rotational strength [10^(-44) esu^2 cm^2] */
  double em_angle;           /* el./mag. dipole trans. moment angle [deg] */
  short group;               /* symmetry group ID */
  };

/* gaussian data struct */
struct gauss_dat {
  char *job_title;                 /* job title */
  char *job_mthd;                  /* calculation method name */
  char *job_basis;                 /* calculation basis set name */
  short job_mthd_id;                /* calculation method ID */
  short job_type_id;                /* calculation type ID */
  int prog_ver;                     /* program version */
  int prog_rev;                     /* program revision */
  long iop[GAUSS_IOP_MAX];          /* array of IOps */
  unsigned n_atoms;                 /* number of atoms */
  unsigned n_geom;                  /* number of structures */
  unsigned n_fdeg;                  /* number of degrees of freedom */
  unsigned n_vib_modes;             /* number of vibrational modes */
  unsigned n_states;                /* number of excited states */
  unsigned n_electrons;             /* number of electrons */
  unsigned n_alpha_electrons;       /* number of alpha electrons */
  unsigned n_beta_electrons;        /* number of beta electrons */
  unsigned n_bg_charges;            /* number of background charges */
  unsigned multiplicity;            /* multiplicity */
  int charge;                       /* charge */
  double **bg_chrg;                 /* background charges */
  struct basis *bs;                 /* basis set */
  struct gauss_atom *atom;          /* atomic data */
  struct gauss_mo *mo_a;            /* alpha MOs */
  struct gauss_mo *mo_b;            /* beta MOs */
  struct gauss_state *state;        /* excited states */
  struct gauss_geom *geom;          /* structures */
  struct gauss_nmr *nmr;            /* nmr shifts */
  struct gauss_freq *freq;          /* vibrational frequencies */
  };

/* data struct for reading basis set from gaussian log file */
struct gauss_log_bs {
  unsigned at_id;                   /* center ID */
  unsigned at_num;                  /* atomic number */
  double coord[3];                  /* center coordinates */
  struct basis_shell *shell;        /* basis set shell data */
  };

/* allocate memory for gaussian data struct */
struct gauss_dat *gauss_dat_new(void);
/* free memory allocated for gaussian data struct */
struct gauss_dat* gauss_dat_free(struct gauss_dat*);

/* allocate memory for gaussian atom data struct */
struct gauss_atom *gauss_atom_new(unsigned);
/* free memory allocated for gaussian atom data struct */
struct gauss_atom* gauss_atom_free(struct gauss_atom*, unsigned);

/* allocate memory for gaussian MO data struct */
struct gauss_mo *gauss_mo_new(unsigned);
/* free memory allocated for gaussian MO data struct */
struct gauss_mo* gauss_mo_free(struct gauss_mo*, unsigned);

/* allocate memory for gaussian state data struct */
struct gauss_state *gauss_state_new(unsigned);
/* free memory allocated for gaussian state data struct */
struct gauss_state* gauss_state_free(struct gauss_state*, unsigned);

/* allocate memory for gaussian geom data struct */
struct gauss_geom *gauss_geom_new(unsigned);
/* free memory allocated for gaussian geom data struct */
struct gauss_geom* gauss_geom_free(struct gauss_geom*, unsigned, unsigned);

/* allocate vector of NMR shift data structs */
struct gauss_nmr *gauss_nmr_new(unsigned);
/* free memory allocated for gaussian nmr data struct */
struct gauss_nmr* gauss_nmr_free(struct gauss_nmr*, unsigned);

/* allocate array of vibrational-frequency data structs */
struct gauss_freq *gauss_freq_new(unsigned);
/* free memory allocated for gaussian frequency data struct array */
struct gauss_freq* gauss_freq_free(struct gauss_freq*, unsigned, unsigned);

/* copy data from one gaussian data struct to another */
void gauss_atom_copy(struct gauss_atom*, struct gauss_atom*);

/* merge data from two gaussian data struct to one */
struct gauss_dat* gauss_dat_merge(struct gauss_dat*, struct gauss_dat*);
/* merge two molecular orbitals data structs */
struct gauss_mo* gauss_mo_merge(struct gauss_mo*, unsigned, unsigned,
  struct gauss_mo*, unsigned, unsigned);

/* convert gaussian data to generic molecular file format */
struct gen_mol* gauss_dat_gen(struct gauss_dat*);
/* convert gaussian data to xyz molecular file format */
struct xyz_mol* gauss_dat_xyz(struct gauss_dat*);

/* return internal code of program version */
int gauss_prg_ver_id(char*);
/* return internal code of program revision */
int gauss_prg_rev_id(char*);

/* convert calculation type symbol to internal ID */
short gauss_job_type_id(char*);
/* convert calculation type internal ID to symbol */
char* gauss_job_type_name(short);

/* analyze control line with computational options and return method ID */
short gauss_mthd_id(char*);

/* set IOp codes from one line specification */
void gauss_iop_set(long*, char*);

/* calculate nuclear repulsion energy */
double gauss_mol_nucrep(struct gauss_dat*);
/* get energy of highest occupied molecular orbital */
double gauss_mol_homo(struct gauss_dat*, unsigned*, short*);
/* get energy of lowest unoccupied molecular orbital */
double gauss_mol_lumo(struct gauss_dat*, unsigned*, short*);
/* return empirical formula of molecular structure */
char *gauss_mol_formula(struct gauss_dat*);
/* print out molecular structure coordinates */
void gauss_mol_coord_print(struct gauss_dat*, char*);
/* print out molecular structure coordinates to file */
void gauss_mol_coord_fprint(struct gauss_dat*, FILE*);
/* print out molecular coordinates of background charges */
void gauss_mol_chrg_print(struct gauss_dat*, char*);
/* print out molecular coordinates of background charges to file */
void gauss_mol_chrg_fprint(struct gauss_dat*, FILE*);
/* print out molecular orbital coefficients */
void gauss_mol_mo_print(struct gauss_dat*, short, short, char*);
/* print out molecular orbital coefficients to file */
void gauss_mol_mo_fprint(struct gauss_dat*, short, short, FILE*);

/* set basis specified basis set for given structure data */
void gauss_bs_set(struct gauss_dat*, struct basis*);
/* set basis basis set center coodinates according to atoms */
void gauss_bs_set_coord(struct gauss_dat*);
/* return basis set shell type ID corresponding to gaussian internal code */
unsigned gauss_bs_shell_type(short);
/* return gaussian internal code of basis set shell type */
short gauss_bs_shell_id(unsigned);

/* allocate memory for gaussian log-file basis-set data structure */
struct gauss_log_bs* gauss_log_bs_new(void);
/* free memory allocate for log-file basis-set data structure */
struct gauss_log_bs* gauss_log_bs_free(struct gauss_log_bs*);

/* read data from gaussian log file */
void gauss_log_read(struct gauss_dat*, char*, short);
/* read molecule specification from gaussian log file */
void gauss_log_read_spec(struct gauss_dat*, char*);
/* read molecule specification from open gaussian log file */
void gauss_log_read_spec_f(struct gauss_dat*, FILE*);
/* read coordinates from gaussian log file */
void gauss_log_read_coord(struct gauss_dat*, short, char*);
/* read coordinates from open gaussian log file */
void gauss_log_read_coord_f(struct gauss_dat*, short, FILE*);
/* read basis set from gaussian log file */
void gauss_log_read_basis(struct gauss_dat*, char*);
/* read basis set from open gaussian log file */
void gauss_log_read_basis_f(struct gauss_dat*, FILE*);
/* read frequencies from gaussian log file */
void gauss_log_read_freq(struct gauss_dat*, char*);
/* read frequencies from open gaussian log file */
void gauss_log_read_freq_f(struct gauss_dat*, FILE*);
/* read NMR shifts from gaussian log file */
void gauss_log_read_nmr(struct gauss_dat*, char*);
/* read NMR shifts from open gaussian log file */
void gauss_log_read_nmr_f(struct gauss_dat*, FILE*);
/* read excited states from gaussian log file */
void gauss_log_read_states(struct gauss_dat*, char*);
/* read excited states from open gaussian log file */
void gauss_log_read_states_f(struct gauss_dat*, FILE*);

/* read data from gaussian formatted checkpoint file */
void gauss_fchk_read(struct gauss_dat*, char*);
/* write gaussian data to formatted checkpoint file */
void gauss_fchk_write(struct gauss_dat*, char*);

/* read data from NBO output files (basis set file + data file) */
void gauss_nbo_read(struct gauss_dat*, char*, char*);

/* read instructions, charge and multiplicity from gaussian input file */
void gauss_inp_read_header(struct gauss_dat*, char*);
/* read instructions, charge and multiplicity from open gaussian input file */
void gauss_inp_read_header_f(struct gauss_dat*, FILE*);
/* read coordinates from gaussian input file */
void gauss_inp_read_coord(struct gauss_dat*, char*);
/* read coordinates from open gaussian input file */
void gauss_inp_read_coord_f(struct gauss_dat*, FILE*);
/* read background charges from gaussian input file */
void gauss_inp_read_chrg(struct gauss_dat*, char*);
/* read background charges from open gaussian input file */
void gauss_inp_read_chrg_f(struct gauss_dat*, FILE*);

/* write vibrational-frequency normal modes to molden-type file */
void gauss_freq_write(struct gauss_dat*, char*);

/* calculate full overlap matrix over molecular orbitals */
void gauss_calc_mat_s(struct gauss_dat*, short, double**);
/* orthogonalize molecular orbitals saved in the gausssian data struct */
void gauss_calc_mo_orth(struct gauss_dat*, short);
/* sort molecular orbitals according to their energies */
void gauss_calc_mo_sort(struct gauss_dat*);

#endif
