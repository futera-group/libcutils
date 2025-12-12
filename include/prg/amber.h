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

#ifndef ZF_LIB_PRG_AMBER_H
#define ZF_LIB_PRG_AMBER_H

#include <stdio.h>

/* symbolic constants */
#define AMBER_TOP_NPOINTER 31      /* number of pointers */
#define AMBER_TOP_CHARGE   18.2223 /* amber charge scale factor */

#define AMBER_BOND_ALL           0 /* all bonds */
#define AMBER_BOND_WITH_H        1 /* bond including hydrogen */
#define AMBER_BOND_WITHOUT_H     2 /* bond without hydrogen */

#define AMBER_POINTER_NATOM      0 /* number of atoms */
#define AMBER_POINTER_NTYPES     1 /* number of atom types */
#define AMBER_POINTER_NBONH      2 /* number of bonds containing H */
#define AMBER_POINTER_NBONA      3 /* number of bonds not containing H */
#define AMBER_POINTER_NTHETH     4 /* number of angles containing H */
#define AMBER_POINTER_NTHETA     5 /* number of angles not containing H */
#define AMBER_POINTER_NPHIH      6 /* number of dihedrals constaining H */
#define AMBER_POINTER_NPHIA      7 /* number of dihedrals not containing H */
#define AMBER_POINTER_NHPARM     8 /* not used */
#define AMBER_POINTER_NPARM      9 /* not used */
#define AMBER_POINTER_NEXT      10 /* number of exluded atoms */
#define AMBER_POINTER_NRES      11 /* number of residui */
#define AMBER_POINTER_NCBONA    12 /* NBONA + no. of constrained bonds */
#define AMBER_POINTER_NCTHETA   13 /* NTHETA + no. of constrained angels */
#define AMBER_POINTER_NCPHITA   14 /* NPHITA + no. of constrained dihedrals */
#define AMBER_POINTER_NBOND     15 /* number of bond types */ 
#define AMBER_POINTER_NANGLE    16 /* number of angle types */
#define AMBER_POINTER_NDIHED    17 /* number of dihedral types */
#define AMBER_POINTER_NATYP     18 /* number of atom types */
#define AMBER_POINTER_NHBOND    19 /* number of 10-12 H-bond types */
#define AMBER_POINTER_PERT      20 /* perturbation info included */
#define AMBER_POINTER_NPBOND    21 /* number of perturbed bonds */
#define AMBER_POINTER_NPANGLE   22 /* number of perturbed angles */
#define AMBER_POINTER_NPDIHED   23 /* number of perturbed dihedrals */
#define AMBER_POINTER_NPGBOND   24 /* number of bonds in perturbed group */
#define AMBER_POINTER_NPGANGLE  25 /* number of angles in perturbed group */
#define AMBER_POINTER_NPGDIHED  26 /* number of dihedrals in perturbed group */
#define AMBER_POINTER_BOX       27 /* solvation box info included */
#define AMBER_POINTER_NMAXRES   28 /* number of atoms in largest residuum */
#define AMBER_POINTER_CAP       29 /* CAP option set */
#define AMBER_POINTER_POL       30 /* polarization info included */

#define AMBER_BOX_LASTSOLUTE     0 /* last solute residuum */
#define AMBER_BOX_NMOL           1 /* number of molecules */
#define AMBER_BOX_FIRSTSOLV      2 /* first solvent molecule */

/* molecular topology struct */
struct amber_top {
  char *version;                /* version of the file format */
  char *title;                  /* title line */
  int *pointers;                /* specificators */
  char **atom_names;            /* FF atom names */
  double *charge;               /* FF partial charges */
  int *atom_num;                /* atomic numbers */
  double *mass;                 /* atomic masses */
  int *atom_type_id;            /* atom types in L-J interactions */
  int *n_excl_atoms;            /* numbers of exluded atoms */
  int *nbond_prm_id;            /* non-bonding parameter indices */
  char **res_names;             /* residuum names */
  int *res_pointer;             /* residuum atom pointers */
  double *bond_frc_const;       /* force constants for bonds */
  double *bond_eq_val;          /* bond equilibrium value */
  double *angle_frc_const;      /* force constants for angles */
  double *angle_eq_val;         /* angle equilibrium value */
  double *dihed_frc_const;      /* force constants for dihedrals */
  double *dihed_period;         /* dihedral periodicities */
  double *dihed_phase;          /* dihedral phases */
  double *solty;                /* not used */
  double *lj_acoeff;            /* Lennard-Jones a coefficient */
  double *lj_bcoeff;            /* Lennard-Jones b coefficient */
  int *bond_h;                  /* atoms in bonds with H */
  int *bond_ah;                 /* atoms in bonds without H */
  int *angle_h;                 /* atoms in angles with H */
  int *angle_ah;                /* atoms in angles without H */
  int *dihed_h;                 /* atoms in dihedrals with H */
  int *dihed_ah;                /* atoms in dihedrals without H */
  double *scale_ee;             /* electrostatic scaling factors */
  double *scale_nb;             /* non-bonding scaling factors */
  int *excl_atom_list;          /* exluded atom list */
  double *hbond_acoeff;         /* 10-12 H-bond a coefficient */
  double *hbond_bcoeff;         /* 10-12 H-bond b coefficient */
  double *hbond_cut;            /* not used */
  char **atom_types;            /* atom types */
  char **tree_chain_class;      /* tree joining types (not used) */
  int *join_array;              /* tree joining indices (not used) */
  int *irotat;                  /* atom rotation (not used) */
  double *radii;                /* VdW radii of atoms */
  double *screen;               /* screening factors */
  double *polar;                /* atomic polarizabilities */
  int *box_pointer;             /* solvent box pointers */
  int *box_natoms;              /* number of atoms per molecule */
  double *box_params;           /* parameters of the solvent box */
  char *radius_set;             /* name of the radius set */
  int pol;                      /* polarizability flag */
  };

/* amber topology pointer struct */
struct amber_pointer {
  int natom;    /* number of atoms */
  int ntype;    /* number of atom types */
  int natyp;    /* number of atom types */
  int nres;     /* number of residui */
  int nbond;    /* number of bond types */
  int nangle;   /* number of angle types */
  int ndihed;   /* number of dihedral types */
  int nbondh;   /* number of bonds containing H */
  int nbonda;   /* number of bonds not containing H */
  int nangleh;  /* number of angles containing H */
  int nanglea;  /* number of angles not containing H */
  int ndihedh;  /* number of dihedrals containing H */
  int ndiheda;  /* number of dihedrals not containing H */
  int next;     /* number of excluded atoms */
  int nhbond;   /* numbef of H-bonds */
  };

/* allocate new topology struct */
struct amber_top* amber_top_new(void);
/* free memory allocated for topology struct */
void amber_top_free(struct amber_top*);

/* initialize topology pointer struct according to pointer array */
void amber_pointer_init(struct amber_pointer*, int*);

/* create bond array with atom indices for given atom */
unsigned *amber_top_bonds_atom(struct amber_top*, unsigned, unsigned*,
  double**, short);
/* create bond array with atom indices for given residuum */
unsigned **amber_top_bonds_res(struct amber_top*, unsigned, unsigned*, 
  double**, short);

/* read topology file */
void amber_top_read(struct amber_top*, char*);
/* write topology file */
void amber_top_write(struct amber_top*, char*);

/* copy one topology struct to another */
struct amber_top *amber_top_copy(struct amber_top*);

/* read one structure from amber trajectory file */
int amber_crd_read_one(double*, unsigned, double*, FILE*);
/* read specified structure from amber trajectory file */
int amber_crd_read_id(double*, unsigned, double*, unsigned, FILE*);
/* read coordinates from amber restart file */
short amber_crd_read_rst(double*, unsigned, double*, double*, char*);

/* wrap atom coordinates into the PBC box, keep compact residui */
void amber_crd_wrap(struct amber_top*, double*, double*);
/* write one structure to the amber trajectory file */
void amber_crd_write_one(double*, unsigned, double*, FILE*);
/* write coordinates to amber restart file */
void amber_crd_write_rst(double*, unsigned, double*, double*, char*);

/* return IDs of first and last+1 atom of specified residuum */
void amber_res_atom_id(struct amber_top*, unsigned, unsigned*, unsigned*);
/* return number of atoms in specified residuum */
unsigned amber_res_natoms(struct amber_top*, unsigned);
/* return ID of residuum for given atom */
unsigned amber_res_id(struct amber_top*, unsigned*, unsigned);
/* check if the given residuum is terminal or not */
short amber_res_is_ter(struct amber_top*, unsigned);
/* return array with IDs of unique residui */
unsigned *amber_res_unique(struct amber_top*, unsigned*);
/* calculate center of mass of specified residuum of the system */
void amber_res_mass_center(struct amber_top*, double*, unsigned, double*);

/* convert amber topology and coordinates to acf molecular file format */
struct acf_mol *amber_top_acf(struct amber_top*, double*);
/* convert amber topology and coordinates to acf molecular file format */
struct acf_mol *amber_top_acf_all(struct amber_top*, double*);

/* convert amber topology and coordinates to apc molecular file format */
struct apc_mol *amber_top_apc(struct amber_top*, double*, unsigned);

/* convert amber topology and coordinates to cif molecular file format */
struct cif_mol *amber_top_cif(struct amber_top*, double*);
/* convert amber topology and coordinates to cif molecular file format */
struct cif_mol *amber_top_cif_r(struct amber_top*, double*, unsigned);
/* convert amber topology and coordinates to cif molecular file format */
struct cif_mol *amber_top_cif_rv(struct amber_top*, double*,
  unsigned*, unsigned);

/* convert amber topology and coordinates to pdb molecular file format */
struct pdb_mol *amber_top_pdb(struct amber_top*, double*);
/* convert amber topology and coordinates to pdb molecular file format */
struct pdb_mol *amber_top_pdb_r(struct amber_top*, double*, unsigned);
/* convert amber topology and coordinates to pdb molecular file format */
struct pdb_mol *amber_top_pdb_rv(struct amber_top*, double*,
  unsigned*, unsigned);

/* convert amber topology and coordinates to xyz molecular file format */
struct xyz_mol *amber_top_xyz(struct amber_top*, double*);
/* convert amber topology and coordinates to xyz molecular file format */
struct xyz_mol *amber_top_xyz_r(struct amber_top*, double*, unsigned);
/* convert amber topology and coordinates to xyz molecular file format */
struct xyz_mol *amber_top_xyz_rv(struct amber_top*, double*,
  unsigned*, unsigned);

#endif
