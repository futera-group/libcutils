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

#include <stdio.h>
#include <string.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* read one line from amber topology file, end program with error if EOF
   
   line - space for read chars  
   what - specify what is reading (appears in error message)
   file - pointer to open file */
char* amber_top_read_line(char *what, FILE *file) {
  char *line;
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of amber topology file while reading %s",
      1,what);
  return(line);
  }

/* locate flag in amber topology file (version 7 format)

   flag - the flag which is searched
   file - pointer to open topology file */
void amber_top_read_flag(char *flag, FILE *file) {
  char *line,sflag[128];
  short i;
  sprintf(sflag,"%%FLAG %s",flag);
  for (i=0; i<2; i++) {
    for (line=str_read_line_new(file); line;
         line=str_free(line),line=str_read_line_new(file)) {
      if (str_sub_bfind(line,sflag)) {
        /* drop the line with format specification */
        line = str_free(line);
        line = amber_top_read_line("format",file);
        line = str_free(line);
        return;
        }
      }
    /* try again from the beginning */
    rewind(file);
    }
  /* flag not found, exit with error */
  msg_error_f("cannot found flag %s in amber topology file",1,flag);
  }

/* read numbers from amber topology file (integer version)

   num  - number of integers which should be read
   flag - format flag in topology file
   what - description what is read
   file - pointer to open topology file */
int *amber_top_read_inums(unsigned num, char *flag, char *what, FILE *file) {
  int *v;
  unsigned i;
  amber_top_read_flag(flag,file);
  v = vec_ialloc(num);
  for (i=0; i<num; i++)
    if (fscanf(file,"%d",&(v[i]))!=1)
      msg_error_f("unexpected error while reading %s in topology file",1,what);
  return(v);
  }

/* read numbers from amber topology file (double version)

   num  - number of integers which should be read
   flag - format flag in topology file
   what - description what is read
   file - pointer to open topology file */
double *amber_top_read_fnums(unsigned num, char *flag, char *what, FILE *file) {
  double *v;
  unsigned i;
  amber_top_read_flag(flag,file);
  v = vec_falloc(num);
  for (i=0; i<num; i++)
    if (fscanf(file,"%lf",&(v[i]))!=1)
      msg_error_f("unexpected error while reading %s in topology file",1,what);
  return(v);
  }

/* read and check version of topology file 

   str  - string where the version will be saved
   file - pointer to open file */
char* amber_top_read_string(char *flag, FILE *file) {
  char *line,*check,*string;
  line = str_read_line_new(file);
  if (!line)
    msg_error("cannot read string from amber topology file",1);
  check = str_sprintf("%%%s",flag);
  if (!str_sub_bfind(line,check))
    msg_error_f("expected flag \"%s\" not found",1,flag);
  str_free(check);
  string = str_copy_new(line+8);
  str_free(line);
  str_rtrim(string);
  return(string);
  }

/* read string from amber topology file 

   str  - the string which is read
   flag - format flag in topology file
   what - description what is read
   file - pointer to open topology file */
char* amber_top_read_str(char *flag, char *what, FILE *file) {
  char *line;
  amber_top_read_flag(flag,file);
  line = amber_top_read_line(what,file);
  str_trim(line);
  return(line);
  }

/* read vector of strings from amber topology file 

   ntot   - number of strings which should be read
   flag  - format flag in topology file
   what  - description what is read
   file  - pointer to open topology file */
char **amber_top_read_vstr(unsigned ntot, char *flag, char *what, FILE *file) {
  unsigned i,j,n,nn,nr,width = 4,nrow = 20;
  char **v,*line;
  amber_top_read_flag(flag,file);
  v = vec_salloc(ntot,width+1);
  /* number of rows in topology file */
  nr = (ntot%nrow==0 ? ntot/nrow : ntot/nrow+1);
  /* read names */
  n = 0;
  for (i=0; i<nr; i++) {
    line = amber_top_read_line(what,file);
    /* number of names on the row */
    nn = (i==nr-1 ? (ntot%nrow==0 ? nrow : ntot%nrow) : nrow);
    for (j=0; j<nn; j++) {
      strncpy(v[n],line+j*width,width);
      v[n][width]='\0';
      n++;
      }
    str_free(line);
    }
  return(v);
  }

/* check if polarizability is included in topology file or not
 
   file - pointer to open amber topology file */
int amber_top_read_pol(FILE *file) {
  char *line;
  int found = 0;
  line = str_ffind_new(file,"POLARIZABILITY");
  found = (line ? 1 : 0);
  str_free(line);
  rewind(file);
  return(found);
  }

/* read topology file
   
   t         - pointer to topology struct 
   file_name - name of the topology file */
void amber_top_read(struct amber_top *t, char *file_name) {
  FILE *file;
  int *d,read_pol = 0;
  struct amber_pointer p;
  /* open the topology file */
  file = file_open(file_name,"r");
  /* check if polarizabilities are in the file */
  read_pol = amber_top_read_pol(file);
  /* header */
  t->version = amber_top_read_string("VERSION",file);
  t->title = amber_top_read_str("TITLE","title",file);
  /* data pointers */
  t->pointers = amber_top_read_inums(AMBER_TOP_NPOINTER,"POINTERS",
    "pointers",file);
  amber_pointer_init(&p,t->pointers);
  /* atoms */
  t->atom_names = 
    amber_top_read_vstr(p.natom,"ATOM_NAME","atom names",file);
  t->charge = 
    amber_top_read_fnums(p.natom,"CHARGE","charges",file);
  t->atom_num = 
    amber_top_read_inums(p.natom,"ATOMIC_NUMBER","atomic numbers",file);
  t->mass = 
    amber_top_read_fnums(p.natom,"MASS","masses",file);
  /* LJ potentials */
  t->atom_type_id = amber_top_read_inums(p.natom,"ATOM_TYPE_INDEX",
    "indices of atom types in L-J interactions",file);
  t->n_excl_atoms = amber_top_read_inums(p.natom,"NUMBER_EXCLUDED_ATOMS",
    "numbers of excluded atoms",file);
  t->nbond_prm_id = amber_top_read_inums(p.ntype*p.ntype,"NONBONDED_PARM_INDEX",
    "non-bonding parameter indices",file);
  /* residui */
  t->res_names = amber_top_read_vstr(p.nres,"RESIDUE_LABEL",
    "residuum names",file);
  t->res_pointer = amber_top_read_inums(p.nres,"RESIDUE_POINTER",
    "residuum atom pointers",file);
  /* bond force constants */
  t->bond_frc_const = amber_top_read_fnums(p.nbond,"BOND_FORCE_CONSTANT",
    "force constants for bonds",file);
  t->bond_eq_val = amber_top_read_fnums(p.nbond,"BOND_EQUIL_VALUE",
    "bond equilibrium values",file);
  /* angle force constants */
  t->angle_frc_const = amber_top_read_fnums(p.nangle,"ANGLE_FORCE_CONSTANT",
    "force constants for angles",file);
  t->angle_eq_val = amber_top_read_fnums(p.nangle,"ANGLE_EQUIL_VALUE",
    "angle equilibrium values",file);
  /* dihedral force constants */
  t->dihed_frc_const = amber_top_read_fnums(p.ndihed,"DIHEDRAL_FORCE_CONSTANT",
    "force constants for dihedrals",file);
  t->dihed_period = amber_top_read_fnums(p.ndihed,"DIHEDRAL_PERIODICITY",
    "dihedral periodicities",file);
  t->dihed_phase = amber_top_read_fnums(p.ndihed,"DIHEDRAL_PHASE",
    "dihedral phases",file);
  /* scaling factors */
  t->scale_ee = amber_top_read_fnums(p.ndihed,"SCEE_SCALE_FACTOR",
    "electrostatic scaling factor",file);
  t->scale_nb = amber_top_read_fnums(p.ndihed,"SCNB_SCALE_FACTOR",
    "non-bonding scaling factor",file);
  t->solty = amber_top_read_fnums(p.natyp,"SOLTY","solty array",file);
  /* JL coefficients */
  t->lj_acoeff = amber_top_read_fnums(p.ntype*(p.ntype+1)/2,
    "LENNARD_JONES_ACOEF","Lennard-Jones a coefficients",file);
  t->lj_bcoeff = amber_top_read_fnums(p.ntype*(p.ntype+1)/2,
    "LENNARD_JONES_BCOEF","Lennard-Jones b coefficients",file);
  /* bonds with / without hydrogens */
  t->bond_h = amber_top_read_inums(3*p.nbondh,"BONDS_INC_HYDROGEN",
    "bonds with H",file);
  t->bond_ah = amber_top_read_inums(3*p.nbonda,"BONDS_WITHOUT_HYDROGEN",
    "bonds without H",file);
  /* angles with / without hydrogens */
  t->angle_h = amber_top_read_inums(4*p.nangleh,"ANGLES_INC_HYDROGEN",
    "angles with H", file);
  t->angle_ah = amber_top_read_inums(4*p.nanglea,"ANGLES_WITHOUT_HYDROGEN",
    "angles without H",file);
  /* dihedrals with / without hydrogens */
  t->dihed_h = amber_top_read_inums(5*p.ndihedh,"DIHEDRALS_INC_HYDROGEN",
    "dihedrals with H",file);
  t->dihed_ah = amber_top_read_inums(5*p.ndiheda,"DIHEDRALS_WITHOUT_HYDROGEN",
    "dihedrals without H",file);
  /* excluded atoms */
  t->excl_atom_list = amber_top_read_inums(p.next,"EXCLUDED_ATOMS_LIST",
    "list of excluded atoms",file);
  /* hydrogen bonds */
  t->hbond_acoeff = amber_top_read_fnums(p.nhbond,"HBOND_ACOEF",
    "\"a\" coefficients of H-bonds",file);
  t->hbond_bcoeff = amber_top_read_fnums(p.nhbond,"HBOND_BCOEF",
    "\"b\" coefficients of H-bonds",file);
  t->hbond_cut = amber_top_read_fnums(p.nhbond,"HBCUT","H-bond cutoffs",file);
  /* graph connectivities */
  t->atom_types = amber_top_read_vstr(p.natom,"AMBER_ATOM_TYPE",
    "atom types",file);
  t->tree_chain_class = amber_top_read_vstr(p.natom,"TREE_CHAIN_CLASSIFICATION",
    "tree joining types",file);
  t->join_array = amber_top_read_inums(p.natom,"JOIN_ARRAY",
    "tree joining indices",file);
  /* rotations */
  t->irotat = amber_top_read_inums(p.natom,"IROTAT",
    "atom rotation indices",file);
  /* solvent */
  if (t->pointers[AMBER_POINTER_BOX]>0) {
    t->box_pointer = amber_top_read_inums(3,"SOLVENT_POINTERS","solvent box"
      " pointers",file);
    t->box_natoms = amber_top_read_inums(t->box_pointer[AMBER_BOX_NMOL],
      "ATOMS_PER_MOLECULE","number of atoms per molecule",file);
    t->box_params = amber_top_read_fnums(4,"BOX_DIMENSIONS","solvent box"
      " parameters",file);
    }
  /* radius set */
  t->radius_set = amber_top_read_str("RADIUS_SET","radius set",file);
  /* capping */
  if (t->pointers[AMBER_POINTER_CAP]>0)
    msg_error("I don't know how to treat file with CAP info, sorry...",1);
  /* perturbed potential */
  if (t->pointers[AMBER_POINTER_PERT]>0) 
    msg_error("I don't know how to treat file with PERT info, sorry...",1);
  if (t->pointers[AMBER_POINTER_POL]>0) 
    msg_error("I don't know how to treat file with POL info, sorry...",1);
  /* atom radii and screening */
  t->radii = amber_top_read_fnums(p.natom,"RADII","VdW radii",file);
  t->screen = amber_top_read_fnums(p.natom,"SCREEN","screening factors",file);
  /* polarization */
  if (read_pol)
    t->polar = amber_top_read_fnums(p.natom,"POLARIZABILITY",
      "polarizabilities",file);
  d = amber_top_read_inums(1,"IPOL","polarizability flag",file);
  t->pol = d[0];
  vec_ifree(d);
  /* close the file */
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
