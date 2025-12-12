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

#ifndef ZF_LIB_PRG_GROMACS_H
#define ZF_LIB_PRG_GROMACS_H

#include <stdio.h>
#include <string.h>
#include <mol/molec.h>

/* -------------------------------------------------------------------------- */

/* particle types */
#define GMX_PTYPE_ATOM   1   /* atom */
#define GMX_PTYPE_SHELL  2   /* shell */
#define GMX_PTYPE_VIRT   3   /* virtual site */

/* bonding terms */
#define GMX_BTERM_BOND   1   /* bonds */
#define GMX_BTERM_ANGL   2   /* angles */
#define GMX_BTERM_DIHE   3   /* proper dihedrals */
#define GMX_BTERM_IMPR   4   /* improper dihedrals */
#define GMX_BTERM_PAIR   5   /* 1-4 pairs */
#define GMX_BTERM_NBPR   6   /* other non-bonded pairs */
#define GMX_BTERM_CNST   7   /* constraints */
#define GMX_BTERM_CMAP   8   /* cmap corrections */
#define GMX_BTERM_STTL   9   /* settles */
#define GMX_BTERM_EXCL  10   /* exclusions */

/* fragment types */
#define GMX_FTYPE_PROT   1   /* protein */
#define GMX_FTYPE_WAT    2   /* water */
#define GMX_FTYPE_ION    3   /* ion */
#define GMX_FTYPE_SURF   4   /* surface */

/* -------------------------------------------------------------------------- */

/* Atomic data structure */
struct gmx_atom {
  short particle;          /* particle type */
  unsigned id;             /* atomic ID */
  unsigned num;            /* atomic number */
  char *name;              /* atom name */
  char *type;              /* atom type */
  double charge;           /* atomic charge */
  double mass;             /* atomic mass */
  double sigma;            /* LJ-potential sigma parameter */
  double epsilon;          /* LJ-potential epsilon parameter */
  };

/* Bonding data structure */
struct gmx_bond {
  short type;              /* functional type */
  unsigned n_atoms;        /* number of atoms */
  unsigned n_parms;        /* number of parameters */
  char **name;             /* atom names */
  unsigned *ia;            /* atom IDs */
  unsigned *ir;            /* residuum IDs */
  double *parm;            /* parameters */
  };

/* Topology data structure */
struct gmx_top {
  unsigned n_atom_types;   /* number of atom types */
  unsigned n_bonds;        /* number of bonds */
  unsigned n_angles;       /* number of angles */
  unsigned n_dihedrals;    /* number of dihedrals */
  unsigned n_impropers;    /* number of improper dihedrals */
  unsigned n_14_pairs;     /* number of 1-4 non-bonding pairs */
  unsigned n_nb_pairs;     /* number of other non-bonding pairs */
  unsigned n_constraints;  /* number of constraints */
  unsigned n_settles;      /* number of settles */
  unsigned n_igb_parms;    /* number of implicit-solvent parameters */
  unsigned n_cmap;         /* number of cmap bondings */
  unsigned n_exclusions;   /* number of excluded exclusions */
  struct gmx_atom *type;   /* array of atom types */
  struct gmx_bond *bond;   /* array of bonds */
  struct gmx_bond *angl;   /* array of angles */
  struct gmx_bond *dihe;   /* array of dihedrals */
  struct gmx_bond *impr;   /* array of improper dihedrals */
  struct gmx_bond *pair;   /* array of 1-4 non-bonding pairs */
  struct gmx_bond *nbpr;   /* array of other non-bonding pairs */
  struct gmx_bond *cnst;   /* array of constraints */
  struct gmx_bond *sttl;   /* array of settles */
  struct gmx_bond *igbp;   /* array of implicit-solvent parameters */
  struct gmx_bond *cmap;   /* array of cmap bonding */
  struct gmx_bond *excl;   /* array of exclusions */
  };

/* Residuum data structure */
struct gmx_res {
  char *name;              /* residuum name */
  unsigned id;             /* residuum ID */
  unsigned n_atoms;        /* number of atoms */
  struct gmx_atom *atom;   /* array of atoms */
  struct gmx_top *top;     /* intra-residual topology terms */
  };

/* Molecular data structure */
struct gmx_mol {
  char *file;              /* name of file containing the data */
  char *name;              /* name of the molecule (ITP keyword) */
  short gen_itp;           /* generate ITP file */
  unsigned n_excl;         /* excluded pairs */
  unsigned n_atoms;        /* number of atoms */
  unsigned n_resids;       /* number of residues */
  struct gmx_res *res;     /* array of residues */
  struct gmx_top *top;     /* inter-residual topology terms */
  };

/* Force-field data structure */
struct gmx_ff {
  char *file;              /* name of file containing the force-field data */
  unsigned n_defs;         /* number of pre-processor define directives */
  char **def;              /* array of pre-processor define directives */
  short nb_func;           /* non-bonded function type */
  short comb_rule;         /* non-bonded combination rule */
  short gen_pairs;         /* generate list of 1-4 pair interactions */
  double scale_lj;         /* Lennard-Jones scaling factor */
  double scale_qq;         /* Coulomb scaling factor */
  struct gmx_top *dat;     /* force-field data */
  };

/* Gromacs structure fragment data */
struct gmx_frag {
  short type;              /* fragment type (protein,water,ion,...) */
  char *name;              /* name of the molecule */
  unsigned n_rep;          /* number of duplicates */
  struct gmx_mol *mol;     /* pointer to molecular data */
  };

/* Gromacs structure and force-field data */
struct gmx_dat {
  char *title;             /* system title */
  unsigned n_atoms;        /* total number of atoms */
  unsigned n_resids;       /* total number of residues */
  unsigned n_mols;         /* number of molecules */
  unsigned n_frags;        /* number of system fragments */
  double **crd;            /* atom coordinates */
  double *box;             /* simulation box */
  struct gmx_ff *ff;       /* force field data */
  struct gmx_mol *mol;     /* array of molecular definions */
  struct gmx_frag *frag;   /* array of system fragments */
  };

/* -------------------------------------------------------------------------- */

/* Allocate memory for array of atoms */
struct gmx_atom* gmx_atom_new(unsigned);
/* Clean memory allocated for atom array */
struct gmx_atom* gmx_atom_free(struct gmx_atom*, unsigned);
/* Clean memory allocated in atom structure */
void gmx_atom_clean(struct gmx_atom*);
/* Copy atomic data from one data structure to another */
void gmx_atom_copy(struct gmx_atom*, struct gmx_atom*);
/* Create copy of atom data structure array */
struct gmx_atom* gmx_atom_copy_new(struct gmx_atom*, unsigned);
/* Compare two atom data structures */
short gmx_atom_compare(struct gmx_atom*, struct gmx_atom*);

/* Allocate memory for array of bonds */
struct gmx_bond* gmx_bond_new(unsigned);
/* Clean memory allocated for bond array */
struct gmx_bond* gmx_bond_free(struct gmx_bond*, unsigned);
/* Clean memory allocated in bond data structure */
void gmx_bond_clean(struct gmx_bond*);
/* Copy data from one bond data structure to another */
void gmx_bond_copy(struct gmx_bond*, struct gmx_bond*);
/* Create copy of bonding data structure array */
struct gmx_bond* gmx_bond_copy_new(struct gmx_bond*, unsigned);

/* Allocate memory for topology data */
struct gmx_top* gmx_top_new(void);
/* Clean memory allocated for topology */
struct gmx_top* gmx_top_free(struct gmx_top*);
/* Clean memory allocated in topology data structure */
void gmx_top_clean(struct gmx_top*);
/* Copy topology data from one data structure to another */
void gmx_top_copy(struct gmx_top*, struct gmx_top*);
/* Create copy of topology data structure */
struct gmx_top* gmx_top_copy_new(struct gmx_top*);

/* Allocate memory for array of residues */
struct gmx_res* gmx_res_new(unsigned);
/* Clean memory allocated for residuum array */
struct gmx_res* gmx_res_free(struct gmx_res*, unsigned);
/* Clean memory allocated in residuum structure */
void gmx_res_clean(struct gmx_res*);
/* Copy residuum data from one data structure to another */
void gmx_res_copy(struct gmx_res*, struct gmx_res*);
/* Create copy of residuum data structure array */
struct gmx_res* gmx_res_copy_new(struct gmx_res*, unsigned);
/* Compare two residuum data structures */
short gmx_res_compare(struct gmx_res*, struct gmx_res*);

/* Allocate memory for array of molecules */
struct gmx_mol* gmx_mol_new(unsigned);
/* Clean memory allocated for molecule array */
struct gmx_mol* gmx_mol_free(struct gmx_mol*, unsigned);
/* Clean memory allocated in molecule structure */
void gmx_mol_clean(struct gmx_mol*);
/* Copy molecule data from one data structure to another */
void gmx_mol_copy(struct gmx_mol*, struct gmx_mol*);
/* Create copy of molecule data structure array */
struct gmx_mol* gmx_mol_copy_new(struct gmx_mol*, unsigned);
/* List searching function comparing molecule file names */
short gmx_mol_cmp_file_wrp_search(void*, void*);
/* List sorting function comparing molecule file names */
short gmx_mol_cmp_file_wrp_sort(void*, void*);
/* Compare two molecular data structures */
short gmx_mol_compare(struct gmx_mol*, struct gmx_mol*);

/* Allocate memory for gromacs force-field data structure */
struct gmx_ff* gmx_ff_new(void);
/* Clean memory allocated for gromacs force-field data structure */
struct gmx_ff* gmx_ff_free(struct gmx_ff*);
/* Clean memory allocated in gromacs force-field data structure */
void gmx_ff_clean(struct gmx_ff*);
/* Copy gromacs data from one force-field data structure to another */
void gmx_ff_copy(struct gmx_ff*, struct gmx_ff*);
/* Create copy of gromacs force-field data structure */
struct gmx_ff* gmx_ff_copy_new(struct gmx_ff*);

/* Allocate memory for array of fragments */
struct gmx_frag* gmx_frag_new(unsigned);
/* Clean memory allocated for fragment array */
struct gmx_frag* gmx_frag_free(struct gmx_frag*, unsigned);
/* Clean memory allocated in fragment structure */
void gmx_frag_clean(struct gmx_frag*);
/* Copy fragment data from one data structure to another */
void gmx_frag_copy(struct gmx_frag*, struct gmx_frag*);
/* Create copy of fragment data structure array */
struct gmx_frag* gmx_frag_copy_new(struct gmx_frag*, unsigned);
/* Compare two fragment data structures */
short gmx_frag_compare(struct gmx_frag*, struct gmx_frag*);

/* Allocate memory for gromacs data structure */
struct gmx_dat* gmx_dat_new(void);
/* Clean memory allocated for gromacs data structure */
struct gmx_dat* gmx_dat_free(struct gmx_dat*);
/* Clean memory allocated in gromacs data structure */
void gmx_dat_clean(struct gmx_dat*);
/* Copy gromacs data from one data structure to another */
void gmx_dat_copy(struct gmx_dat*, struct gmx_dat*);
/* Create copy of gromacs data structure */
struct gmx_dat* gmx_dat_copy_new(struct gmx_dat*);

/* -------------------------------------------------------------------------- */

/* Call pre-processor on given file and send data to open pipe */
FILE *gmx_read_pipe_open(char*, char**, unsigned);

/* Add new terms to bonding data array */
void gmx_bond_add(struct gmx_bond**, unsigned*, struct gmx_bond*, unsigned);
/* Delete specified term from bonding data array */
void gmx_bond_del_one(struct gmx_bond**, unsigned, unsigned*);
/* Delete specified terms from bonding data array */
void gmx_bond_del(struct gmx_bond**, unsigned*, unsigned, unsigned*);

/* Return keyword name of topology bonding group */
char *gmx_bond_key_name(short);
/* Return brief description of topology bonding group */
char *gmx_bond_key_desc(short);
/* Return number of atoms involved in definition of specified bonding term */
unsigned gmx_bond_key_natoms(short);

/* Add new atoms to residuum data structure */
void gmx_res_atom_add(struct gmx_res*, struct gmx_atom*, unsigned);
/* Return ID of specified atom in the residuum */
unsigned gmx_res_atom_get_id(struct gmx_res*, char*);

/* Return pointer to specified bonding term */
struct gmx_bond** gmx_top_bond_get(struct gmx_top*, short, unsigned**);
/* Remove all bonding terms involving specified atom */
void gmx_top_bond_del(struct gmx_top*, unsigned, unsigned, short);
/* Sort save dihedral types / definition to proper / improper */
void gmx_top_dihe_prop_impr(struct gmx_top*);

/* Return internal ID of particle type */
short gmx_ff_particle_type_id(char*);
/* Convert internal ID to particle-type symbol */
char *gmx_ff_particle_type_name(short);

/* Add new atoms to specific residuum in molecular structure */
void gmx_mol_atom_add(struct gmx_mol*, unsigned, struct gmx_atom*, unsigned);
/* Return residuum ID of specified atom in a molecule */
unsigned gmx_mol_res_get_id(struct gmx_mol*, unsigned*, unsigned);
/* Set residuum IDs and sort bonding terms to inter/intra residual ones */
void gmx_mol_set_bonding(struct gmx_mol*);

/* Read one molecule definition from open gromacs topology file */
char* gmx_mol_fread_itp_one(struct gmx_mol*, char*, FILE*);
/* Read one molecule definition from gromacs topology file */
void gmx_mol_read_itp_one(struct gmx_mol*, char*, char**, unsigned);
/* Read all molecule definitions from open gromacs topology file */
struct gmx_mol* gmx_mol_fread_itp(unsigned*, FILE*);
/* Read all molecule definitions from gromacs topology file */
struct gmx_mol* gmx_mol_read_itp(char*, unsigned*, char**, unsigned);
/* Write molecular type to open gromacs topology file */
void gmx_mol_fwrite_itp(struct gmx_mol*, double*, FILE*);

/* Update number of atoms and residues in the system by fragment definitions */
void gmx_dat_update_nums(struct gmx_dat*);

/* Return system ID of specified atom in the system */
unsigned gmx_dat_atom_get_id(struct gmx_dat*, unsigned, unsigned, 
  unsigned, unsigned);

/* Calculate COM of given residuum */
void gmx_dat_res_com(struct gmx_dat*, struct gmx_res*, unsigned, double*);
/* Return IDs of residui within given distance of specific residui */
unsigned* gmx_dat_res_get_list_dr(struct gmx_dat*, unsigned*, unsigned,
  double, unsigned*);

/* Add molecular topologies to the system */
void gmx_dat_mol_add(struct gmx_dat*, struct gmx_mol*, unsigned);
/* Return pointer to molecular data structure with specified name */
struct gmx_mol* gmx_dat_mol_get(struct gmx_dat*, char*);

/* Add molecular coordinates to the system */
void gmx_dat_frag_add(struct gmx_dat*, struct gmx_frag*, unsigned,
  double**, unsigned);
/* Delete specified molecular fragments from the system */
void gmx_dat_frag_del(struct gmx_dat*, unsigned*, unsigned);
/* Return fragment ID of specified atom in the system */
unsigned gmx_dat_frag_get_id(struct gmx_dat*, unsigned*, unsigned*, unsigned);
/* Compact fragment array in gromacs data structure (use repetitions) */
void gmx_dat_frag_compact(struct gmx_dat*);
/* Expand fragment array in gromacs data structure (no repetitions) */
void gmx_dat_frag_expand(struct gmx_dat*);

/* Read structure from open G96 file */
int gmx_dat_fread_g96(struct gmx_dat*, long unsigned*, double*, short*,
  short, FILE*);
/* Read structure from G96 file */
void gmx_dat_read_g96(struct gmx_dat*, char*);
/* Read structure from open GRO file */
void gmx_dat_fread_gro(struct gmx_dat*, FILE*);
/* Read structure from GRO file */
void gmx_dat_read_gro(struct gmx_dat*, char*);
/* Read structure from open XYZ file */
void gmx_dat_fread_xyz(struct gmx_dat*, FILE*);
/* Read structure from XYZ file */
void gmx_dat_read_xyz(struct gmx_dat*, char*);
/* Read structure from open PDB file */
void gmx_dat_fread_pdb(struct gmx_dat*, FILE*);
/* Read structure from PDB file */
void gmx_dat_read_pdb(struct gmx_dat*, char*);
/* Read data from gromacs topology file */
void gmx_dat_read_top(struct gmx_dat*, char*);
/* Read data-block keyword from line extracted from Gromacs file */
char *gmx_dat_read_top_key(char*);

/* Write structure to open G96 file */
void gmx_dat_fwrite_g96(struct gmx_dat*, long unsigned, double, short, FILE*);
/* Write structure to G96 file */
void gmx_dat_write_g96(struct gmx_dat*, char*);
/* Write structure to open GRO file */
void gmx_dat_fwrite_gro(struct gmx_dat*, FILE*);
/* Write structure to GRO file */
void gmx_dat_write_gro(struct gmx_dat*, char*);
/* Write structure to open gromacs topology file */
void gmx_dat_fwrite_top(struct gmx_dat*, short, double*, FILE*);
/* Write structure to gromacs topology file */
void gmx_dat_write_top(struct gmx_dat*, short, double*, char*);

/* Convert gromacs data structure to XYZ molecular format */
struct xyz_mol* gmx_dat_conv_to_xyz(struct gmx_dat*);
/* Convert gromacs data structure to XYZ molecular format */
struct xyz_mol* gmx_dat_conv_to_xyz_rv(struct gmx_dat*,
  unsigned*, unsigned);
/* Convert gromacs data structure to PDB molecular format */
struct pdb_mol* gmx_dat_conv_to_pdb(struct gmx_dat*);
/* Convert gromacs data structure to PDB molecular format */
struct pdb_mol* gmx_dat_conv_to_pdb_rv(struct gmx_dat*,
  unsigned*, unsigned);

#endif
