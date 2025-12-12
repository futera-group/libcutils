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

#ifndef ZF_LIB_PRG_DLPOLY_H
#define ZF_LIB_PRG_DLPOLY_H

/* units */
#define	DLPOLY_UNIT_EV   1    /* electron volts */
#define	DLPOLY_UNIT_KC   2    /* k-calories per mol */
#define	DLPOLY_UNIT_KJ   3    /* k-joules per mol */
#define	DLPOLY_UNIT_KB   4    /* kelvin per boltzmann */
#define	DLPOLY_UNIT_IN   5    /* internal units (10 k-joules per mol) */

/* tethering potentials */
#define DLPOLY_TETH_HM   1    /* harmonic */	
#define DLPOLY_TETH_RH   2    /* restraint */
#define DLPOLY_TETH_QR   3    /* quartic */

/* bond potentials */
#define DLPOLY_BOND_HM   1    /* harmonic */
#define DLPOLY_BOND_MR   2    /* morse */
#define DLPOLY_BOND_12   3    /* 12-6 */
#define DLPOLY_BOND_LJ   4    /* lennard-jones */
#define DLPOLY_BOND_RH   5    /* restraint */
#define DLPOLY_BOND_QR   6    /* quartic */
#define DLPOLY_BOND_BC   7    /* buckingham */
#define DLPOLY_BOND_CL   8    /* coulomb */
#define DLPOLY_BOND_SF   9    /* shifted FENE */
#define DLPOLY_BOND_AM  10    /* AMOEBA */

/* angle potentials */
#define DLPOLY_ANGL_HM   1    /* harmonic */
#define DLPOLY_ANGL_QR   2    /* quartic */
#define DLPOLY_ANGL_TH   3    /* truncated harmonic */
#define DLPOLY_ANGL_SH   4    /* screened harmonic */
#define DLPOLY_ANGL_SV   5    /* screened vessal */
#define DLPOLY_ANGL_TV   6    /* truncated vessal */
#define DLPOLY_ANGL_HC   7    /* harmonic cosine */
#define DLPOLY_ANGL_CS   8    /* cosine */
#define DLPOLY_ANGL_MM   9    /* MM3 stretch-bend */
#define DLPOLY_ANGL_SS  10    /* compass stretch-stretch */
#define DLPOLY_ANGL_SB  11    /* compass stretch-bend */
#define DLPOLY_ANGL_CM  12    /* compass all terms */
#define DLPOLY_ANGL_AM  13    /* AMOEBA */
#define DLPOLY_ANGL_KK  14    /* KKY */

/* dihedral potentials */
#define DLPOLY_DIHE_CS   1    /* cosine */
#define DLPOLY_DIHE_HM   2    /* harmonic */
#define DLPOLY_DIHE_HC   3    /* harmonic cosine */
#define DLPOLY_DIHE_TC   4    /* triple cosine */
#define DLPOLY_DIHE_RB   5    /* ryckaert-bellemans */
#define DLPOLY_DIHE_FR   6    /* fluorinated ryckaert-bellemans */
#define DLPOLY_DIHE_OP   7    /* OPLS torsion */

/* inversion potentials */
#define DLPOLY_INVR_HM   1    /* harmonic */
#define DLPOLY_INVR_HC   2    /* harmonic cosine */
#define DLPOLY_INVR_PL   3    /* planar */
#define DLPOLY_INVR_EP   4    /* extended planar */
#define DLPOLY_INVR_CL   5    /* calcite */

/* pair potentials */
#define DLPOLY_PAIR_TB   1    /* tabulation */
#define DLPOLY_PAIR_12   2    /* 12-6 */
#define DLPOLY_PAIR_LJ   3    /* lennard-jones */
#define DLPOLY_PAIR_NM   4    /* n-m */
#define DLPOLY_PAIR_BC   5    /* buckingham */
#define DLPOLY_PAIR_BH   6    /* born-huggins-meyer */
#define DLPOLY_PAIR_HB   7    /* 12-10 h-bond */
#define DLPOLY_PAIR_SF   8    /* shifted force n-m */
#define DLPOLY_PAIR_MR   9    /* morse */
#define DLPOLY_PAIR_WC  10    /* shifted weeks-chandler-anderson */
#define DLPOLY_PAIR_DP  11    /* standard DPD (groot-warren) */
#define DLPOLY_PAIR_AM  12    /* AMOEBA 14-7 */

/* metal potentials */
#define DLPOLY_METL_EA   1    /* EAM */
#define DLPOLY_METL_EE   2    /* EEAM */
#define DLPOLY_METL_2B   3    /* 2BEAM */
#define DLPOLY_METL_2E   4    /* 2BEEAM */
#define DLPOLY_METL_FS   5    /* Finnis-Sinclair */
#define DLPOLY_METL_EF   6    /* Extended Finnis-Sinclair */
#define DLPOLY_METL_SC   7    /* Sutton-Chen */
#define DLPOLY_METL_GP   8    /* Gupta */
#define DLPOLY_METL_MB   9    /* MBPC */

/* tersoff potentials */
#define DLPOLY_TERS_TR   1    /* TERS-type potential */
#define DLPOLY_TERS_KH   2    /* KIHS-type potential */

/* three-body potentials */
#define DLPOLY_TBPT_HM   1    /* harmonic */
#define DLPOLY_TBPT_TH   2    /* truncated harmonic */
#define DLPOLY_TBPT_SH   3    /* screened harmonic */
#define DLPOLY_TBPT_SV   4    /* screened vessal */
#define DLPOLY_TBPT_TV   5    /* truncated vessal */
#define DLPOLY_TBPT_HB   6    /* h-bonds */

/* four-body potentials */
#define DLPOLY_FBPT_HM   1    /* harmonic */
#define DLPOLY_FBPT_HC   2    /* harmonic cosine */
#define DLPOLY_FBPT_PL   3    /* planar */

/* external field */
#define DLPOLY_EFLD_EL   1    /* electric field */
#define DLPOLY_EFLD_OS   2    /* oscillating shear */
#define DLPOLY_EFLD_CS   3    /* continuous shear */
#define DLPOLY_EFLD_GV   4    /* gravitational field */
#define DLPOLY_EFLD_MG   5    /* magnetic field */
#define DLPOLY_EFLD_SP   6    /* containing sphere */
#define DLPOLY_EFLD_RW   7    /* repulsive wall */
#define DLPOLY_EFLD_XP   8    /* x-piston */
#define DLPOLY_EFLD_HA   9    /* molecule in HR zone */
#define DLPOLY_EFLD_HB  10    /* HR zone (pull out) */
#define DLPOLY_EFLD_HC  11    /* HR zone (pull in) */
#define DLPOLY_EFLD_OE  12    /* oscillating electric field */

/* DL_POLY atomic data struct */
struct dlpoly_atom {
  char *name;                 /* name of the atom */
  unsigned n_atoms;           /* number of this atoms in the molecule */
  short frozen;               /* frozen / free atoms */
  double mass;                /* atomic mass */
  double charge;              /* atomic charge */
  };

/* DL_POLY force-field shell parameters */
struct dlpoly_prm_shell {
  unsigned id1;               /* site index of core */
  unsigned id2;               /* site index of shell */
  double k2;                  /* force constant of core-shell spring */
  double k4;                  /* quartic (anharmonic) force constant */
  };

/* DL_POLY force-field constrain parameters */
struct dlpoly_prm_cnst {
  unsigned id1;               /* first atomic site index */
  unsigned id2;               /* second atomic site index */
  double bl;                  /* constraint bond length */
  };

/* DL_POLY force-field PMF bond-length parameters */
struct dlpoly_prm_pmf {
  unsigned id1;               /* first atomic site index */
  unsigned id2;               /* second atomic site index */
  double w1;                  /* first atomic site weight */
  double w2;                  /* second atomic site weight */
  };

/* DL_POLY force-field rigid parameters */
struct dlpoly_prm_rig {
  unsigned n_sites;           /* number of sites in rigid unit */
  unsigned *site;             /* site atomic IDs */
  };

/* DL_POLY force-field thethering parameters */
struct dlpoly_prm_teth {
  short pot;                  /* potential ID */
  unsigned site;              /* atomic site ID */
  double v1;                  /* potential parameter #1 */
  double v2;                  /* potential parameter #2 */
  };

/* DL_POLY force-field bond parameters */
struct dlpoly_prm_bond {
  short pot;                  /* potential ID */
  unsigned id1;               /* first atomic site ID */
  unsigned id2;               /* second atomic site ID */
  double v1;                  /* potential parameter #1 */
  double v2;                  /* potential parameter #2 */
  double v3;                  /* potential parameter #3 */
  double v4;                  /* potential parameter #4 */
  };

/* DL_POLY force-field angle parameters */
struct dlpoly_prm_angle {
  short pot;                  /* potential ID */
  unsigned id1;               /* first atomic site ID */
  unsigned id2;               /* second atomic site ID */
  unsigned id3;               /* third atomic site ID */
  double v1;                  /* potential parameter #1 */
  double v2;                  /* potential parameter #2 */
  double v3;                  /* potential parameter #3 */
  double v4;                  /* potential parameter #4 */
  };

/* DL_POLY force-field dihedral parameters */
struct dlpoly_prm_dihed {
  short pot;                  /* potential ID */
  unsigned id1;               /* first atomic site ID */
  unsigned id2;               /* second atomic site ID */
  unsigned id3;               /* third atomic site ID */
  unsigned id4;               /* fouth atomic site ID */
  double v1;                  /* potential parameter #1 */
  double v2;                  /* potential parameter #2 */
  double v3;                  /* potential parameter #3 */
  double v4;                  /* potential parameter #4 */
  double v5;                  /* potential parameter #5 */
  double v6;                  /* potential parameter #6 */
  double v7;                  /* potential parameter #7 */
  };

/* DL_POLY force-field inversion parameters */
struct dlpoly_prm_inv {
  short pot;                  /* potential ID */
  unsigned id1;               /* first atomic site ID */
  unsigned id2;               /* second atomic site ID */
  unsigned id3;               /* third atomic site ID */
  unsigned id4;               /* fouth atomic site ID */
  double v1;                  /* potential parameter #1 */
  double v2;                  /* potential parameter #2 */
  double v3;                  /* potential parameter #3 */
  };

/* DL_POLY molecular data struct */
struct dlpoly_mol {
  char *name;                     /* name of the molecule */
  unsigned n_molecules;           /* number of this molecules in the system */
  unsigned n_atoms;               /* number of atoms */
  struct dlpoly_atom *atom;       /* array of atoms */
  unsigned n_shells;              /* number of shells */
  struct dlpoly_prm_shell *shell; /* array of shells */
  unsigned n_constraints;         /* number of constraints */
  struct dlpoly_prm_cnst *cnst;   /* array of constraints */
  unsigned n_pmfs;                /* number of PMF bond-lengths */
  struct dlpoly_prm_pmf *pmf;     /* array of PMF bond-lengths */
  unsigned n_rigids;              /* number of rigid units */
  struct dlpoly_prm_rig *rig;     /* array of rigid units */
  unsigned n_teths;               /* number of tethered atoms */
  struct dlpoly_prm_teth *teth;   /* array of tethered atoms */
  unsigned n_bonds;               /* number of bonds */
  struct dlpoly_prm_bond *bond;   /* array of bonds */
  unsigned n_angles;              /* number of angles */
  struct dlpoly_prm_angle *angle; /* array of angles */
  unsigned n_dihedrals;           /* number of dihedral angles */
  struct dlpoly_prm_dihed *dihed; /* array of dihedral angles */
  unsigned n_inversions;          /* number of inversions */
  struct dlpoly_prm_inv *inv;     /* array of inversions */
  };

/* DL_POLY force-field van der Waals parameters */
struct dlpoly_prm_vdw {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  short pot;                      /* potential ID */
  double v1;                      /* potential parameter #1 */
  double v2;                      /* potential parameter #2 */
  double v3;                      /* potential parameter #3 */
  double v4;                      /* potential parameter #4 */
  double v5;                      /* potential parameter #5 */
  };

/* DL_POLY force-field metal parameters */
struct dlpoly_prm_metal {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  short pot;                      /* potential ID */
  double v1;                      /* potential parameter #1 */
  double v2;                      /* potential parameter #2 */
  double v3;                      /* potential parameter #3 */
  double v4;                      /* potential parameter #4 */
  double v5;                      /* potential parameter #5 */
  double v6;                      /* potential parameter #6 */
  double v7;                      /* potential parameter #7 */
  double v8;                      /* potential parameter #8 */
  double v9;                      /* potential parameter #9 */
  };

/* DL_POLY force-field RDF parameters */
struct dlpoly_prm_rdf {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  };

/* DL_POLY force-field tersoff parameters */
struct dlpoly_prm_ters {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  short pot;                      /* potential ID */
  short type;                     /* term type (self of cross) */
  double *v;                      /* potential parameters */
  };

/* DL_POLY force-field three-body-potential parameters */
struct dlpoly_prm_tbp {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  char *at3;                      /* third atom type */
  short pot;                      /* potential ID */
  double v1;                      /* potential parameter #1 */
  double v2;                      /* potential parameter #2 */
  double v3;                      /* potential parameter #3 */
  double v4;                      /* potential parameter #4 */
  double v5;                      /* potential parameter #5 */
  };

/* DL_POLY force-field four-body-potential parameters */
struct dlpoly_prm_fbp {
  char *at1;                      /* first atom type */
  char *at2;                      /* second atom type */
  char *at3;                      /* third atom type */
  char *at4;                      /* four atom type */
  short pot;                      /* potential ID */
  double v1;                      /* potential parameter #1 */
  double v2;                      /* potential parameter #2 */
  double v3;                      /* potential parameter #3 */
  };

/* DL_POLY force-field external-field parameters */
struct dlpoly_prm_efld {
  short pot;                      /* potential ID */
  double v1;                      /* potential parameter #1 */
  double v2;                      /* potential parameter #2 */
  double v3;                      /* potential parameter #3 */
  double v4;                      /* potential parameter #4 */
  double v5;                      /* potential parameter #5 */
  };

/* DL_POLY force-field data struct */
struct dlpoly_fld {
  char *header;                   /* header line */
  short units;                    /* FF parameter units */
  unsigned n_molecules;           /* number of molecules */
  struct dlpoly_mol *mol;         /* array of molecules */
  unsigned n_vdws;                /* number of VDW potentials */
  struct dlpoly_prm_vdw *vdw;     /* array of VDW potentials */
  unsigned n_metals;              /* number of metal potentials */
  struct dlpoly_prm_metal *metal; /* array of metal potentials */
  unsigned n_rdfs;                /* number of RDF pairs */
  struct dlpoly_prm_rdf *rdf;     /* array of RDF pairs */
  unsigned n_tersoffs;            /* number of tersoff potentials */
  struct dlpoly_prm_ters *ters;   /* array of tersoff potentials */
  unsigned n_tbps;                /* number of three-body potentials */
  struct dlpoly_prm_tbp *tbp;     /* array of three-body potentials */
  unsigned n_fbps;                /* number of four-body potentials */
  struct dlpoly_prm_fbp *fbp;     /* array of four-body potentials */
  unsigned n_ext_fields;          /* number of external-field potentials */
  struct dlpoly_prm_efld *efld;   /* external field */
  };

/* DL_POLY structure data */
struct dlpoly_cfg {
  char *header;                   /* header line */
  short data;                     /* data-type ID */
  short pbc;                      /* periodic-bounday ID */
  unsigned n_atoms;               /* number of atoms */
  double **cell;                  /* cell vectors */
  char **sym;                     /* atom symbols */
  double **crd;                   /* coordinates */
  double **vel;                   /* velocities */
  double **frc;                   /* forces */
  };

/* allocate memory for DL_POLY atomic data */
struct dlpoly_atom *dlpoly_atom_new(unsigned);
/* free memory allocated for DL_POLY atomic data */
struct dlpoly_atom *dlpoly_atom_free(struct dlpoly_atom*, unsigned);

/* allocate memory for DL_POLY molecular data */
struct dlpoly_mol *dlpoly_mol_new(unsigned);
/* free memory allocated for DL_POLY molecular data */
struct dlpoly_mol *dlpoly_mol_free(struct dlpoly_mol*, unsigned);

/* allocate memory for DL_POLY angle parameters */
struct dlpoly_prm_angle *dlpoly_prm_angle_new(unsigned);
/* free memory allocated for DL_POLY angle-potential data */
struct dlpoly_prm_angle *dlpoly_prm_angle_free(struct dlpoly_prm_angle*);
/* allocate memory for DL_POLY bond parameters */
struct dlpoly_prm_bond *dlpoly_prm_bond_new(unsigned);
/* free memory allocated for DL_POLY bond-potential data */
struct dlpoly_prm_bond *dlpoly_prm_bond_free(struct dlpoly_prm_bond*);
/* allocate memory for DL_POLY constrain parameters */
struct dlpoly_prm_cnst *dlpoly_prm_cnst_new(unsigned);
/* free memory allocated for DL_POLY constrain-parameter data */
struct dlpoly_prm_cnst *dlpoly_prm_cnst_free(struct dlpoly_prm_cnst*);
/* allocate memory for DL_POLY dihedral angle parameters */
struct dlpoly_prm_dihed *dlpoly_prm_dihed_new(unsigned);
/* free memory allocated for DL_POLY dihedral angle-potential data */
struct dlpoly_prm_dihed *dlpoly_prm_dihed_free(struct dlpoly_prm_dihed*);
/* allocate memory for DL_POLY external-field parameters */
struct dlpoly_prm_efld *dlpoly_prm_efld_new(unsigned);
/* free memory allocated for DL_POLY exteranal-field potential data */
struct dlpoly_prm_efld *dlpoly_prm_efld_free(struct dlpoly_prm_efld*);
/* allocate memory for DL_POLY four-body parameters */
struct dlpoly_prm_fbp *dlpoly_prm_fbp_new(unsigned);
/* free memory allocated for DL_POLY four-body potential data */
struct dlpoly_prm_fbp *dlpoly_prm_fbp_free(struct dlpoly_prm_fbp*, unsigned);
/* allocate memory for DL_POLY inversion parameters */
struct dlpoly_prm_inv *dlpoly_prm_inv_new(unsigned);
/* free memory allocated for DL_POLY inversion-potential data */
struct dlpoly_prm_inv *dlpoly_prm_inv_free(struct dlpoly_prm_inv*);
/* allocate memory for DL_POLY metal parameters */
struct dlpoly_prm_metal *dlpoly_prm_metal_new(unsigned);
/* free memory allocated for DL_POLY metal potential data */
struct dlpoly_prm_metal *dlpoly_prm_metal_free(struct dlpoly_prm_metal*,
   unsigned);
/* allocate memory for DL_POLY PMF parameters */
struct dlpoly_prm_pmf *dlpoly_prm_pmf_new(unsigned n);
/* free memory allocated for DL_POLY PMF data */
struct dlpoly_prm_pmf *dlpoly_prm_pmf_free(struct dlpoly_prm_pmf*);
/* allocate memory for DL_POLY RDF parameters */
struct dlpoly_prm_rdf *dlpoly_prm_rdf_new(unsigned);
/* free memory allocated for DL_POLY RDF pair data */
struct dlpoly_prm_rdf *dlpoly_prm_rdf_free(struct dlpoly_prm_rdf*, unsigned);
/* allocate memory for DL_POLY rigid parameters */
struct dlpoly_prm_rig *dlpoly_prm_rig_new(unsigned);
/* free memory allocated for DL_POLY rigid-unit data */
struct dlpoly_prm_rig *dlpoly_prm_rig_free(struct dlpoly_prm_rig*, unsigned);
/* allocate memory for DL_POLY shell parameters */
struct dlpoly_prm_shell *dlpoly_prm_shell_new(unsigned);
/* free memory allocated for DL_POLY shell-parameter data */
struct dlpoly_prm_shell *dlpoly_prm_shell_free(struct dlpoly_prm_shell*);
/* allocate memory for DL_POLY three-body parameters */
struct dlpoly_prm_tbp *dlpoly_prm_tbp_new(unsigned);
/* free memory allocated for DL_POLY three-body potential data */
struct dlpoly_prm_tbp *dlpoly_prm_tbp_free(struct dlpoly_prm_tbp*, unsigned);
/* allocate memory for DL_POLY tersoff parameters */
struct dlpoly_prm_ters *dlpoly_prm_ters_new(unsigned);
/* free memory allocated for DL_POLY tersoff potential data */
struct dlpoly_prm_ters *dlpoly_prm_ters_free(struct dlpoly_prm_ters*, unsigned);
/* allocate memory for DL_POLY tethering parameters */
struct dlpoly_prm_teth *dlpoly_prm_teth_new(unsigned);
/* free memory allocated for DL_POLY tethering-potential data */
struct dlpoly_prm_teth *dlpoly_prm_teth_free(struct dlpoly_prm_teth*);
/* allocate memory for DL_POLY VdW parameters */
struct dlpoly_prm_vdw *dlpoly_prm_vdw_new(unsigned);
/* free memory allocated for DL_POLY VdW potential data */
struct dlpoly_prm_vdw *dlpoly_prm_vdw_free(struct dlpoly_prm_vdw*, unsigned);

/* allocate memory for DL_POLY force-field parameters */
struct dlpoly_fld *dlpoly_fld_new(void);
/* free memory allocated for DL_POLY force-field potential data  */
struct dlpoly_fld *dlpoly_fld_free(struct dlpoly_fld*);
/* return internal ID of force-field units */
short dlpoly_fld_unit_id(char*);
/* convert internal ID of force-field units to its symbol */
char *dlpoly_fld_unit_sym(short);
/* return number of atoms defined in the force field */
unsigned dlpoly_fld_atom_n(struct dlpoly_fld*);
/* find specific atom in the force-field data */
short dlpoly_fld_atom_find(struct dlpoly_fld*, char*, unsigned*, unsigned*);
/* read force-field parameters from the external file */
void dlpoly_fld_read(struct dlpoly_fld*, char*);
/* write force-field parameters to external file */
void dlpoly_fld_write(struct dlpoly_fld*, char*);

/* return internal ID of angle potential */
short dlpoly_prm_angle_pot_id(char*);
/* return number of variables of specific angle potential */
unsigned dlpoly_prm_angle_pot_nvar(short);
/* convert internal ID of angle potential to its symbol */
char *dlpoly_prm_angle_pot_sym(short, short);
/* return internal ID of bond potential */
short dlpoly_prm_bond_pot_id(char*);
/* return number of variables of specific bond potential */
unsigned dlpoly_prm_bond_pot_nvar(short);
/* convert internal ID of bond potential to its symbol */
char *dlpoly_prm_bond_pot_sym(short, short);
/* return internal ID of dihedral angle potential */
short dlpoly_prm_dihed_pot_id(char*);
/* return number of variables of specific dihedral angle potential */
unsigned dlpoly_prm_dihed_pot_nvar(short);
/* convert internal ID of dihedral angle potential to its symbol */
char *dlpoly_prm_dihed_pot_sym(short);
/* return internal ID of external-field potential */
short dlpoly_prm_efld_pot_id(char*);
/* convert internal ID of external-field potential to its symbol */
char *dlpoly_prm_efld_pot_sym(short);
/* return internal ID of four-body potential */
short dlpoly_prm_fbp_pot_id(char*);
/* convert internal ID of four-body potential to its symbol */
char *dlpoly_prm_fbp_pot_sym(short);
/* return internal ID of inversion potential */
short dlpoly_prm_inv_pot_id(char*);
/* convert internal ID of inversion potential to its symbol */
char *dlpoly_prm_inv_pot_sym(short);
/* return internal ID of metal potential */
short dlpoly_prm_metal_pot_id(char*);
/* convert internal ID of metal potential to its symbol */
char *dlpoly_prm_metal_pot_sym(short);
/* return internal ID of three-body potential */
short dlpoly_prm_tbp_pot_id(char*);
/* convert internal ID of three-body potential to its symbol */
char *dlpoly_prm_tbp_pot_sym(short);
/* return internal ID of tersoff potential */
short dlpoly_prm_ters_pot_id(char*);
/* convert internal ID of tersoff potential to its symbol */
char *dlpoly_prm_ters_pot_sym(short);
/* return internal ID of tethering potential */
short dlpoly_prm_teth_pot_id(char*);
/* convert internal ID of tethering potential to its symbol */
char *dlpoly_prm_teth_pot_sym(short);
/* return internal ID of VdW potential */
short dlpoly_prm_vdw_pot_id(char*);
/* return number of variables of VDW potential */
unsigned dlpoly_prm_vdw_pot_nvar(short);
/* convert internal ID of VdW potential to its symbol */
char *dlpoly_prm_vdw_pot_sym(short);

/* allocate memory for DL_POLY structure data */
struct dlpoly_cfg *dlpoly_cfg_new(void);
/* free memory allocated for DL_POLY structure data */
struct dlpoly_cfg *dlpoly_cfg_free(struct dlpoly_cfg*);
/* read structure data from the external file */
void dlpoly_cfg_read(struct dlpoly_cfg*, char*);
/* write structure data to external file */
void dlpoly_cfg_write(struct dlpoly_cfg*, char*);

/* convert DL_POLY structure & force-field data to XYZ molecule format */
struct xyz_mol *dlpoly_to_xyz(struct dlpoly_cfg*, struct dlpoly_fld*);

#endif
