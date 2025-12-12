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

#ifndef ZF_LIB_PRG_LAMMPS_H
#define ZF_LIB_PRG_LAMMPS_H

#include <stdio.h>
#include "prg/dlpoly.h"

/* atomic types / styles */
#define LAMMPS_ATOM_ANGL   1  /* bead-spring polymers with stiffness */
#define LAMMPS_ATOM_ATOM   2  /* coarse-grain liquids, solids, metals */
#define LAMMPS_ATOM_BODY   3  /* arbitrary bodies */
#define LAMMPS_ATOM_BOND   4  /* bead-spring polymers */
#define LAMMPS_ATOM_CHRG   5  /* atomic system with charges */
#define LAMMPS_ATOM_DIPL   6  /* system with dipolar particles */
#define LAMMPS_ATOM_DPDP   7  /* DPD particles */
#define LAMMPS_ATOM_ELEC   8  /* electronic force field */
#define LAMMPS_ATOM_ELLP   9  /* aspherical particles */
#define LAMMPS_ATOM_FULL  10  /* bio-molecules */
#define LAMMPS_ATOM_LINE  11  /* rigid bodies */
#define LAMMPS_ATOM_MESO  12  /* SPH particles */
#define LAMMPS_ATOM_MOLE  13  /* uncharged molecules */
#define LAMMPS_ATOM_PERI  14  /* mesoscopic peridynamic models */ 
#define LAMMPS_ATOM_SMDP  15  /* solid and fluid SPH particles */
#define LAMMPS_ATOM_SPHR  16  /* granular models */
#define LAMMPS_ATOM_TEMP  17  /* small molecules with fixed topology */
#define LAMMPS_ATOM_TRIP  18  /* rigid bodies */
#define LAMMPS_ATOM_WPCK  19  /* AWPMD */
#define LAMMPS_ATOM_HYBR  20  /* hybrid */

/* potential types */
#define LAMMPS_POT_PAIR    1  /* pair */
#define LAMMPS_POT_PRIJ    2  /* pair IJ */
#define LAMMPS_POT_BOND    3  /* bond */
#define LAMMPS_POT_ANGL    4  /* angle */
#define LAMMPS_POT_DIHE    5  /* dihedral */
#define LAMMPS_POT_IMPR    6  /* improper */
#define LAMMPS_POT_BNBN    7  /* bond-bond */
#define LAMMPS_POT_BNAN    8  /* bond-angle */
#define LAMMPS_POT_MBTR    9  /* middle bond-torsion */
#define LAMMPS_POT_EBTR   10  /* end bond-torsion */
#define LAMMPS_POT_ANTR   11  /* angle-torsion */
#define LAMMPS_POT_ANAN   12  /* angle-angle */
#define LAMMPS_POT_AATR   13  /* angle-angle-torsion */
#define LAMMPS_POT_BB13   14  /* bond-bond 13 */

/* bond types */
#define LAMMPS_BOND_NONE   1  /* turn off bonded iteractions */
#define LAMMPS_BOND_HYBR   2  /* multiple bond types */
#define LAMMPS_BOND_CLS2   3  /* COMPASS 2 */
#define LAMMPS_BOND_FENE   4  /* FENE */
#define LAMMPS_BOND_FEEX   5  /* FENE with variable size particles */
#define LAMMPS_BOND_HARM   6  /* harmonic */
#define LAMMPS_BOND_MORS   7  /* Morse */
#define LAMMPS_BOND_NONL   8  /* nonlinear */
#define LAMMPS_BOND_QUAR   9  /* breakable quartic */
#define LAMMPS_BOND_TABL  10  /* tabulated */

/* angle types */
#define LAMMPS_ANGL_NONE   1  /* turn off angle iteractions */
#define LAMMPS_ANGL_ZERO   2  /* topology without interactions */
#define LAMMPS_ANGL_HYBR   3  /* multiple angle types */
#define LAMMPS_ANGL_CHRM   4  /* CHARMM */
#define LAMMPS_ANGL_CLS2   5  /* COMPASS 2 */
#define LAMMPS_ANGL_COSN   6  /* cosine of angle */
#define LAMMPS_ANGL_COSD   7  /* difference of angle cosines */
#define LAMMPS_ANGL_COSP   8  /* DREIDING angle cosine */
#define LAMMPS_ANGL_COSS   9  /* squared angle cosine */
#define LAMMPS_ANGL_HARM  10  /* harmonic */
#define LAMMPS_ANGL_TABL  11  /* tabulated */

/* diheral-angle types */
#define LAMMPS_DIHE_NONE   1  /* turn off dihedral interactions */
#define LAMMPS_DIHE_ZERO   2  /* topology without interactions */
#define LAMMPS_DIHE_HYBR   3  /* multiple dihedral-angle types */
#define LAMMPS_DIHE_CHRM   4  /* CHARMM */
#define LAMMPS_DIHE_CLS2   5  /* COMPASS 2 */
#define LAMMPS_DIHE_HARM   6  /* harmonic */
#define LAMMPS_DIHE_HELX   7  /* helix */
#define LAMMPS_DIHE_MULT   8  /* multi-harmonic */
#define LAMMPS_DIHE_OPLS   9  /* OPLS */

/* pair-potential types */
#define LAMMPS_PAIR_NONE   1  /* turn off pairwise interactions */
#define LAMMPS_PAIR_ZERO   2  /* neighbor list without interactions */
#define LAMMPS_PAIR_HYBR   3  /* multiple style pairwise interaction */
#define LAMMPS_PAIR_HBOV   4  /* multiple style of superposed pairwise int. */
#define LAMMPS_PAIR_ADPT   5  /* angular dependent potential */
#define LAMMPS_PAIR_AIRB   6  /* AIREBO with LJ */
#define LAMMPS_PAIR_AIRM   7  /* AIREBO with Morse */
#define LAMMPS_PAIR_BECK   8  /* Beck potential */
#define LAMMPS_PAIR_BODY   9  /* interaction between body particles */
#define LAMMPS_PAIR_BOPT  10  /* BOP potential */
#define LAMMPS_PAIR_BORN  11  /* Born-Mayer-Huggins potential */
#define LAMMPS_PAIR_BRCL  12  /* Born-Mayer-Huggins with long-range Coulomb */
#define LAMMPS_PAIR_BRCC  13  /* Born-Mayer-Huggins with Coulomb and cores */
#define LAMMPS_PAIR_BRCM  14  /* Born-Mayer-Huggins with MSM */
#define LAMMPS_PAIR_BRCW  15  /* Born-Mayer-Huggins with Wolf potential */
#define LAMMPS_PAIR_BRWN  16  /* Brownian potential */
#define LAMMPS_PAIR_BRPL  17  /* Brownian potential with polydispersity */
#define LAMMPS_PAIR_BUCK  18  /* Buckingham potential */
#define LAMMPS_PAIR_BCCT  19  /* Buckingham potential with cutoff Coulomb */
#define LAMMPS_PAIR_BCCL  20  /* Buckingham potential with long-range Coulomb */
#define LAMMPS_PAIR_BCCC  21  /* Buckingham potential with Coulomb and core */
#define LAMMPS_PAIR_BCCM  22  /* Buckingham potential with MSM */
#define LAMMPS_PAIR_BCLC  23  /* Buckingham potential with long-range Coulomb */
#define LAMMPS_PAIR_CLLD  24  /* Integrated colloidal potential */
#define LAMMPS_PAIR_COMB  25  /* charge-optimized many body (COMB) */
#define LAMMPS_PAIR_CMB3  26  /* COMB3 */
#define LAMMPS_PAIR_CLCT  27  /* cutoff Coulomb */
#define LAMMPS_PAIR_CLDB  28  /* cutoff Coulomb with Debye screening */
#define LAMMPS_PAIR_CLDS  29  /* Coulomb via damped shifted forces */
#define LAMMPS_PAIR_CLLG  30  /* long-range Coulomb */
#define LAMMPS_PAIR_CLLC  31  /* long-range Coulomb with core */
#define LAMMPS_PAIR_CLMS  32  /* long-range MSM Coulomb */
#define LAMMPS_PAIR_CLSM  33  /* Coulombics via Streitz/Mintmire Slater orb. */
#define LAMMPS_PAIR_CLWF  34  /* Coulombics via Wolf potential */
#define LAMMPS_PAIR_DPDP  35  /* dissipative particle dynamics (DPD) */
#define LAMMPS_PAIR_DPDT  36  /* DPD thermostatting */
#define LAMMPS_PAIR_DSMC  37  /* direct simulation monte carlo (DSMC) */
#define LAMMPS_PAIR_EAMP  38  /* embedded atom method (EAM) */
#define LAMMPS_PAIR_EAMA  39  /* alloy EAM */
#define LAMMPS_PAIR_EAMF  40  /* Finnis-Sinclair EAM */
#define LAMMPS_PAIR_EIMP  41  /* embedde ion method (EIM) */
#define LAMMPS_PAIR_GAUS  42  /* Gaussian potential */
#define LAMMPS_PAIR_GBEP  43  /* Gay-Berne ellipsoidal potential */
#define LAMMPS_PAIR_GRHR  44  /* granular potential with Hertzian interaction */
#define LAMMPS_PAIR_GRHK  45  /* granular potential with history effects */
#define LAMMPS_PAIR_GRHH  46  /* granular potential without history effects */
#define LAMMPS_PAIR_HBDL  47  /* DREIDING hydrogen bonding LJ potential */
#define LAMMPS_PAIR_HBDM  48  /* DREIDING hydrogen bonding Morse potential */
#define LAMMPS_PAIR_KIMP  49  /* KIM project interface */
#define LAMMPS_PAIR_CBOP  50  /* long-range bond-order potential (LCBOP) */
#define LAMMPS_PAIR_LJPT  51  /* LJ potential */
#define LAMMPS_PAIR_CHCC  52  /* CHARMM potential with cutoff Coulomb */
#define LAMMPS_PAIR_CHCI  53  /* CHARMM potential for implicit solvent */
#define LAMMPS_PAIR_CHCL  54  /* CHARMM potential with long-range Coulomb */
#define LAMMPS_PAIR_CHCM  55  /* CHARMM potential with long-range MSM Coulomb */
#define LAMMPS_PAIR_CLS2  56  /* COMPASS2 with no Coulomb */
#define LAMMPS_PAIR_C2CC  57  /* COMPASS2 with cutoff Coulomb */
#define LAMMPS_PAIR_C2LC  58  /* COMPASS2 with long-range Coulomb */
#define LAMMPS_PAIR_LJCB  59  /* LJ with cubic after inflection point */
#define LAMMPS_PAIR_LJCT  60  /* LJ with cutoff and no Coulomb */
#define LAMMPS_PAIR_LJCC  61  /* LJ with cutoff Coulomb */
#define LAMMPS_PAIR_LJCD  62  /* LJ with Debye screening added to Coulomb */
#define LAMMPS_PAIR_LJCS  63  /* LJ with damped-shifted-forces Coulomb */
#define LAMMPS_PAIR_LJCL  64  /* LJ with long-range Coulomb */
#define LAMMPS_PAIR_LJCR  65  /* LJ with long-range Coulomb and core */
#define LAMMPS_PAIR_LJCM  66  /* LJ with long-range MSM Coulomb */
#define LAMMPS_PAIR_LJDC  67  /* point dipoles with cutoff */
#define LAMMPS_PAIR_LJDL  68  /* point dipoles with long-range Ewald */
#define LAMMPS_PAIR_LJT4  69  /* LJ with cutoff Coulomb for TIP4P water */
#define LAMMPS_PAIR_LJL4  70  /* LJ with long-range Coulomb for TIP4P water */
#define LAMMPS_PAIR_LJVS  71  /* LJ for variable size particles */
#define LAMMPS_PAIR_GXLJ  72  /* GROMACS LJ */
#define LAMMPS_PAIR_GXCL  73  /* GROMACS LJ and Coulomb */
#define LAMMPS_PAIR_LLCL  74  /* long-range LJ and long-range Coulomb */
#define LAMMPS_PAIR_LLDL  75  /* long-range LJ and long-range point dipoles */
#define LAMMPS_PAIR_LLT4  76  /* long-range LJ and long-range TIP4P */
#define LAMMPS_PAIR_SMLJ  77  /* smoothed LJ potential */
#define LAMMPS_PAIR_SMLN  78  /* linear smoothed LJ potential */
#define LAMMPS_PAIR_LJ96  79  /* LJ 9/6 potential */
#define LAMMPS_PAIR_LBPT  80  /* hydrodynamic lubrication forces */
#define LAMMPS_PAIR_LBPL  81  /* hydr. lubrication forces with polydispesity */
#define LAMMPS_PAIR_LBFD  82  /* fast lubrication dynamics */
#define LAMMPS_PAIR_LBFP  83  /* fast lubrication dynamics with polydisp. */
#define LAMMPS_PAIR_MEAM  84  /* modified embeddded atom method (MEAM) */
#define LAMMPS_PAIR_MIEP  85  /* Mie potential */
#define LAMMPS_PAIR_MORS  86  /* Morse potential */
#define LAMMPS_PAIR_NB3B  87  /* nonbonded 3-body harmonic potential */
#define LAMMPS_PAIR_NMPT  88  /* NM potential */
#define LAMMPS_PAIR_NMCC  89  /* NM potential with cutoff Coulomb */
#define LAMMPS_PAIR_NMCL  90  /* NM potential with long-range Coulomb */
#define LAMMPS_PAIR_PREP  91  /* peridynamic EPS potential */
#define LAMMPS_PAIR_PRLP  92  /* peridynamic LPS potential */
#define LAMMPS_PAIR_PRPM  93  /* peridynamic PMB potential */
#define LAMMPS_PAIR_PRVS  94  /* peridynamic VES potential */
#define LAMMPS_PAIR_PLMR  95  /* polymorphic 3-body potential */
#define LAMMPS_PAIR_REAX  96  /* ReaxFF potential */
#define LAMMPS_PAIR_REBO  97  /* 2nd generation REBO potential of Brenner */
#define LAMMPS_PAIR_RESQ  98  /* Everaers RE-Squared ellipsoidal potential */
#define LAMMPS_PAIR_SNAP  99  /* SNAP quantum-accurate potential */
#define LAMMPS_PAIR_SOFT 100  /* Soft cosine potential */
#define LAMMPS_PAIR_SWPT 101  /* Stillinger-Weber 3-body potential */
#define LAMMPS_PAIR_TABL 102  /* tabulated pair potential */
#define LAMMPS_PAIR_TERS 103  /* Tersoff 3-body potential */
#define LAMMPS_PAIR_TERM 104  /* modified Tersoff 3-body potential */
#define LAMMPS_PAIR_TERZ 105  /* Tersoff/ZBL 3-body potential */
#define LAMMPS_PAIR_TP4C 106  /* Coulomb for TIP4P water with/withoug LJ */
#define LAMMPS_PAIR_TP4L 107  /* long-range Coulomb for TIP4P water */
#define LAMMPS_PAIR_TRIL 108  /* LJ potential between triangles */
#define LAMMPS_PAIR_VASH 109  /* Vashishta 2-body and 3-body potential */
#define LAMMPS_PAIR_YUKW 110  /* Yukawa potential */
#define LAMMPS_PAIR_YUKC 111  /* screened Yukawa potential */
#define LAMMPS_PAIR_ZBLP 112  /* Ziegler-Biersack-Littmark potential */

/* LAMMPS potential data */
struct lammps_pot {
  unsigned type;                  /* potential type ID */
  short style;                    /* potential style (harmonic,quartic,...) */
  unsigned n_ids;                 /* number of atom IDs */
  unsigned *id;                   /* array of atom IDs */
  unsigned n_coeffs;              /* number of coefficients */
  double *coeff;                  /* array of coefficients */
  };

/* LAMMPS particle type data */
struct lammps_type {
  char *name;                     /* atomic symbol */
  unsigned num;                   /* atomic number */
  double mass;                    /* mass */
  };

/* LAMMPS atomic data */
struct lammps_atom {
  unsigned type;                  /* atomic type ID */
  unsigned mol;                   /* molecular ID */
  double charge;                  /* atomic charge */
  double *crd;                    /* coordinates */
  double *vel;                    /* velocities */
  double *rvl;                    /* radial velocity (electrons) */
  double *avl;                    /* angular velocity (spheres) */
  double *ang;                    /* angular momemtum (ellipsoids) */
  };

/* LAMMPS simulation box data */
struct lammps_box {
  short *pbc;                     /* periodicity indicator (x,y,z) */
  double *min;                    /* lower boundaries (x,y,z) */
  double *max;                    /* upper boundaries (x,y,z) */
  double *tilt;                   /* titlting factors (xy,xz,yz) */
  double **vector;                /* cell vectors */
  };

/* LAMMPS trajectory data */
struct lammps_trj {
  char title[2][81];              /* title and remark lines */
  short has_box;                  /* trajectory contains box info */
  short has_4d;                   /* trajectory contains 4D data */
  short has_chrg;                 /* trajectory contains charges */
  long unsigned n_frames;         /* number of frames in the trajectory */
  long unsigned n_steps;          /* number of MD steps */
  unsigned save_freq;             /* trajecotry saving frequency */
  double timestep;                /* MD integration time step */
  };

/* LAMMPS structure and force-field data */
struct lammps_dat {
  char *header;                   /* title line */
  unsigned n_atoms;               /* number of atoms */
  unsigned n_bonds;               /* number of bonds */
  unsigned n_angles;              /* number of valence angles */
  unsigned n_dihedrals;           /* number of dihedral angles */
  unsigned n_impropers;           /* number of improper angles */
  unsigned n_atom_types;          /* number of atomic types */
  unsigned n_bond_types;          /* number of bond types */
  unsigned n_angle_types;         /* number of valence-angle types */
  unsigned n_dihed_types;         /* number of dihedral-angle types */
  unsigned n_improp_types;        /* number of improper-angle types */
  unsigned n_extra_bonds;         /* number of extra bonds */
  unsigned n_extra_angles;        /* number of extra valence angles */
  unsigned n_extra_dihedrals;     /* number of extra dihedral angles */
  unsigned n_extra_impropers;     /* number of extra improper angles */
  unsigned n_extra_specials;      /* number of extra special sites */
  unsigned n_ellipsoids;          /* number of defined ellipsoids */
  unsigned n_lines;               /* number of defined lines */
  unsigned n_triangles;           /* number of defined triangles */
  unsigned n_bodies;              /* number of defined bodies */
  short atom_style;               /* atom specification style */
  struct lammps_box *box;         /* simulation box */
  struct lammps_atom *atom;       /* array of atoms */
  struct lammps_pot *bond;        /* array of interatomic bonds */
  struct lammps_pot *angle;       /* array of valence angles */
  struct lammps_pot *dihed;       /* array of dihedral angles */
  struct lammps_pot *impr;        /* array of imporper angles */
  struct lammps_type *type;       /* array of particle types */
  struct lammps_pot *pot_bond;    /* bond potentials */
  struct lammps_pot *pot_angle;   /* angle potentials */
  struct lammps_pot *pot_dihed;   /* dihedral-angle potentials */
  struct lammps_pot *pot_pair;    /* pair potentials */
  struct lammps_trj *trj;         /* trajectory data */
  };

/* allocate memory for LAMMPS atomic data */
struct lammps_atom *lammps_atom_new(unsigned);
/* free memory allocated for LAMMPS atomic data */
struct lammps_atom *lammps_atom_free(struct lammps_atom*, unsigned);
/* return internal ID of atom style */
short lammps_atom_type_id(char*);
/* return symbol of atom specification type */
char* lammps_atom_type_sym(short);
/* return number of variables / coefficients in specific atom style */
unsigned lammps_atom_type_nvar(short);
/* return number of velocity coefficients in specific atom style */
unsigned lammps_atom_vel_nvar(short);
/* read list of atoms from LAMMPS data file */
struct lammps_atom* lammps_atom_read(short, unsigned, unsigned, FILE*);
/* read atomic velocities from LAMMPS data file */
void lammps_atom_read_vel(struct lammps_atom*, short, unsigned,
  unsigned, FILE*);
/* write list of atoms to LAMMPS data file */
void lammps_atom_write(struct lammps_atom*, short, unsigned, FILE*);
/* write list of atomic velocities to LAMMPS data file */
void lammps_atom_write_vel(struct lammps_atom*, short, unsigned, FILE*);

/* allocate memory for LAMMPS simulation box data */
struct lammps_box *lammps_box_new(void);
/* free memory allocated for LAMMPS simulation box data */
struct lammps_box *lammps_box_free(struct lammps_box*);
/* set box parameters from given side lengths and angle cosines */
void lammps_box_set(struct lammps_box*, double, double, double,
  double, double, double);
/* calculate cell vectors from cell side lengts and tilting factors */
void lammps_box_set_vectors(struct lammps_box*);
/* calculate simulation-box side lengths */
void lammps_box_get_side(struct lammps_box*, double*, double*, double*);
/* calculate simulation-box side lengths (vector format) */
void lammps_box_get_side_v(struct lammps_box*, double*);
/* calculate simulation-box angles */
void lammps_box_get_angle(struct lammps_box*, double*, double*, double*);
/* calculate simulation-box angles (vector format) */
void lammps_box_get_angle_v(struct lammps_box*, double*);
/* write simulation box parameters to LAMMPS data file */
void lammps_box_write(struct lammps_box*, FILE*);

/* allocate memory for LAMMPS structure and force-field data */
struct lammps_dat *lammps_dat_new(void);
/* free memory allocated for LAMMPS structure and force-field data */
struct lammps_dat *lammps_dat_free(struct lammps_dat*);
/* read structure and force-field data from LAMMPS data file */
void lammps_dat_read(struct lammps_dat*, char*);
/* write structure and force-field data to LAMMPS data file */
void lammps_dat_write(struct lammps_dat*, char*);

/* allocate memory for LAMMPS particle type data */
struct lammps_type *lammps_type_new(unsigned);
/* free memory allocated for LAMMPS particle type data  */
struct lammps_type *lammps_type_free(struct lammps_type*, unsigned);
/* read atomic types from LAMMPS data file */
struct lammps_type *lammps_type_read(unsigned, FILE*);
/* write atomic types to LAMMPS data file */
void lammps_type_write(struct lammps_type*, unsigned, FILE*);

/* allocate memory for LAMMPS potential data */
struct lammps_pot *lammps_pot_new(unsigned);
/* free memory allocated for LAMMPS potential data */
struct lammps_pot *lammps_pot_free(struct lammps_pot*, unsigned);
/* return internal ID of potential type */
short lammps_pot_type_id(char*);
/* read potential specification from LAMMPS data file */
struct lammps_pot* lammps_pot_read(short, short, short, short, unsigned, 
  unsigned, unsigned, FILE*);
/* write potential specification to LAMMPS data file */
void lammps_pot_write(struct lammps_pot*, short, short, unsigned, FILE*);

/* return number of variables / coefficients in specific bond potential */
unsigned lammps_pot_bond_nvar(short);
/* return internal ID of bond potential */
short lammps_pot_bond_id(char*);
/* return symbol of bond potential */
char* lammps_pot_bond_sym(short);

/* return number of variables / coefficients in specific angle potential */
unsigned lammps_pot_angle_nvar(short);
/* return internal ID of angle potential */
short lammps_pot_angle_id(char*);
/* return symbol of angle potential */
char* lammps_pot_angle_sym(short);

/* return number of variables / coefficients in specific dihedral potential */
unsigned lammps_pot_dihed_nvar(short);
/* return internal ID of dihedral-angle potential */
short lammps_pot_dihed_id(char*);
/* return symbol of dihedral-angle potential */
char* lammps_pot_dihed_sym(short);

/* return internal ID of pair potential */
short lammps_pot_pair_id(char*);
/* return symbol of pair potential */
char* lammps_pot_pair_sym(short);

/* convert DLPOLY data to LAMMPS format */
struct lammps_dat *lammps_from_dlpoly(struct dlpoly_cfg*, struct dlpoly_fld*);
/* convert XYZ molecular data to LAMMPS format */
struct lammps_dat* lammps_from_xyz(struct xyz_mol*);
/* convert LAMMPS data to CIF crystallographic format */
struct cif_mol *lammps_to_cif(struct lammps_dat*);
/* convert LAMMPS data to XYZ molecular format */
struct xyz_mol *lammps_to_xyz(struct lammps_dat*);

/* allocate memory for LAMMPS trajectory data */
struct lammps_trj *lammps_trj_new(void);
/* free memory allocated for LAMMPS trajectory data */
struct lammps_trj *lammps_trj_free(struct lammps_trj*);
/* read header of the Lammps DCD trajectory file */
FILE* lammps_trj_dcd_read_header(struct lammps_dat*, char*);
/* read one structure from Lammps DCD trajectory file */
int lammps_trj_dcd_read_frame(struct lammps_dat*, float*, FILE*);

#endif
