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
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

#define FLAG_I 0
#define FLAG_F 1
#define FLAG_S 2

/* -------------------------------------------------------------------------- */

/* write flag and format to topology file

   flag - the data flag
   frmt - the format
   file - open file stream */
void amber_top_write_flag(char *flag, char *frmt, FILE *file) {
  char line[82];
  unsigned i;
  /* flag */
  for (i=0; i<80; i++)
    line[i] = ' ';
  line[0] = '%';
  line[1] = 'F';
  line[2] = 'L';
  line[3] = 'A';
  line[4] = 'G';
  for (i=0; (i+6)<80 && i<str_length(flag); i++)
    line[i+6] = flag[i];
  line[80] = '\n';
  line[81] = '\0';
  fprintf(file,"%s",line);
  /* format */
  for (i=0; i<80; i++)
    line[i] = ' ';
  line[0] = '%';
  line[1] = 'F';
  line[2] = 'O';
  line[3] = 'R';
  line[4] = 'M';
  line[5] = 'A';
  line[6] = 'T';
  line[7] = '(';
  for (i=0; (i+8)<80 && i<str_length(frmt); i++)
    line[i+8] = frmt[i];
  if ((i+8)<80)
    line[i+8] = ')';
  line[80] = '\n';
  line[81] = '\0';
  fprintf(file,"%s",line);
  }

/* write version info to topology file

   version - the version string
   file    - open file stream */
void amber_top_write_version(char *version, FILE *file) {
  char line[82];
  unsigned i;
  for (i=0; i<80; i++)
    line[i] = ' ';
  line[0] = '%';
  line[1] = 'V';
  line[2] = 'E';
  line[3] = 'R';
  line[4] = 'S';
  line[5] = 'I';
  line[6] = 'O';
  line[7] = 'N';
  for (i=0; (i+8)<80 && i<str_length(version); i++)
    line[i+8] = version[i];
  line[80] = '\n';
  line[81] = '\0';
  fprintf(file,"%s",line);
  }

/* write single-value data to topology file

   val  - the data value
   flag - amber topology flag
   type - specify type of data
   file - pointer to open topology file */
void amber_top_write_val(void *val, char *flag, short type, FILE *file) {
  char frmt[10];
  switch (type) {
    case FLAG_I: sprintf(frmt,"1I8"); break;
    case FLAG_F: sprintf(frmt,"1E16.8"); break;
    case FLAG_S: sprintf(frmt,"1a80"); break;
    }
  amber_top_write_flag(flag,frmt,file);
  switch (type) {
    case FLAG_I: fprintf(file,"%8d\n",*((int*)val)); break;
    case FLAG_F: fprintf(file,"%16.8E\n",*((double*)val)); break;
    case FLAG_S: fprintf(file,"%-80.80s\n",(char*)val); break;
    }
  }

/* write data from vector to topology file

   vec  - the vector including data
   num  - number of values in the vector
   flag - amber topology flag
   type - specify type of data
   file - pointer to open topology file */
void amber_top_write_vec(void *vec, unsigned num, char *flag, 
  short type, FILE *file) {
  char frmt[10];
  switch (type) {
    case FLAG_I: sprintf(frmt,"10I8"); break;
    case FLAG_F: sprintf(frmt,"5E16.8"); break;
    case FLAG_S: sprintf(frmt,"20a4"); break;
    }
  amber_top_write_flag(flag,frmt,file);
  if (num>0) {
    switch (type) {
      case FLAG_I: vec_ifprint(file,(int*)vec,num,10,8); break;
      case FLAG_F: vec_ffprint(file,(double*)vec,num,5,16,8,VEC_FRM_E); break;
      case FLAG_S: vec_sfprint(file,(char**)vec,num,20,4); break;
      }
    }
  else
    fprintf(file,"\n");
  }

/* write topology file

  t         - pointer to topology struct
  file_name - name of the topology file */
void amber_top_write(struct amber_top *t, char *file_name) {
  FILE *file;
  struct amber_pointer p;
  /* output stream */
  if (file_name)
    file = file_open(file_name,"w");
  else
    file = stdout;
  /* data pointers */
  amber_pointer_init(&p,t->pointers);
  /* header */
  amber_top_write_version(t->version,file);
  amber_top_write_flag("TITLE","20a4",file);
  fprintf(file,"%-80.80s\n",t->title);
  /* data pointers */
  amber_top_write_vec(t->pointers,AMBER_TOP_NPOINTER,"POINTERS",FLAG_I,file);
  /* atoms */
  amber_top_write_vec(t->atom_names,p.natom,"ATOM_NAME",FLAG_S,file);
  amber_top_write_vec(t->charge,p.natom,"CHARGE",FLAG_F,file);
  amber_top_write_vec(t->atom_num,p.natom,"ATOMIC_NUMBER",FLAG_I,file);
  amber_top_write_vec(t->mass,p.natom,"MASS",FLAG_F,file);
  amber_top_write_vec(t->atom_type_id,p.natom,"ATOM_TYPE_INDEX",FLAG_I,file);
  amber_top_write_vec(t->n_excl_atoms,p.natom,"NUMBER_EXCLUDED_ATOMS",
    FLAG_I,file);
  /* LJ parameters */
  amber_top_write_vec(t->nbond_prm_id,p.ntype*p.ntype,"NONBONDED_PARM_INDEX",
    FLAG_I,file);
  /* residui */
  amber_top_write_vec(t->res_names,p.nres,"RESIDUE_LABEL",FLAG_S,file);
  amber_top_write_vec(t->res_pointer,p.nres,"RESIDUE_POINTER",FLAG_I,file);
  /* bond force constants */
  amber_top_write_vec(t->bond_frc_const,p.nbond,"BOND_FORCE_CONSTANT",
    FLAG_F,file);
  amber_top_write_vec(t->bond_eq_val,p.nbond,"BOND_EQUIL_VALUE",FLAG_F,file);
  /* angle force constants */
  amber_top_write_vec(t->angle_frc_const,p.nangle,"ANGLE_FORCE_CONSTANT",
    FLAG_F,file);
  amber_top_write_vec(t->angle_eq_val,p.nangle,"ANGLE_EQUIL_VALUE",
    FLAG_F,file);
  /* dihedral force constants */
  amber_top_write_vec(t->dihed_frc_const,p.ndihed,"DIHEDRAL_FORCE_CONSTANT",
    FLAG_F,file);
  amber_top_write_vec(t->dihed_period,p.ndihed,"DIHEDRAL_PERIODICITY",
    FLAG_F,file);
  amber_top_write_vec(t->dihed_phase,p.ndihed,"DIHEDRAL_PHASE",FLAG_F,file);
  /* scaling factors */
  amber_top_write_vec(t->scale_ee,p.ndihed,"SCEE_SCALE_FACTOR",FLAG_F,file);
  amber_top_write_vec(t->scale_nb,p.ndihed,"SCNB_SCALE_FACTOR",FLAG_F,file);
  amber_top_write_vec(t->solty,p.natyp,"SOLTY",FLAG_F,file);
  /* LJ parameters */
  amber_top_write_vec(t->lj_acoeff,p.ntype*(p.ntype+1)/2,"LENNARD_JONES_ACOEF",
    FLAG_F,file);
  amber_top_write_vec(t->lj_bcoeff,p.ntype*(p.ntype+1)/2,"LENNARD_JONES_BCOEF",
    FLAG_F,file);
  /* bonds with / without hydrogens */
  amber_top_write_vec(t->bond_h,3*p.nbondh,"BONDS_INC_HYDROGEN",FLAG_I,file);
  amber_top_write_vec(t->bond_ah,3*p.nbonda,"BONDS_WITHOUT_HYDROGEN",
    FLAG_I,file);
  /* angles with / without hydrogens */
  amber_top_write_vec(t->angle_h,4*p.nangleh,"ANGLES_INC_HYDROGEN",FLAG_I,file);
  amber_top_write_vec(t->angle_ah,4*p.nanglea,"ANGLES_WITHOUT_HYDROGEN",
    FLAG_I,file);
  /* dihedrals with / without hydrogens */
  amber_top_write_vec(t->dihed_h,5*p.ndihedh,"DIHEDRALS_INC_HYDROGEN",
    FLAG_I,file);
  amber_top_write_vec(t->dihed_ah,5*p.ndiheda,"DIHEDRALS_WITHOUT_HYDROGEN",
    FLAG_I,file);
  /* exluded atoms */
  amber_top_write_vec(t->excl_atom_list,p.next,"EXCLUDED_ATOMS_LIST",
    FLAG_I,file);
  /* hydrogen bonds */
  amber_top_write_vec(t->hbond_acoeff,p.nhbond,"HBOND_ACOEF",FLAG_F,file);
  amber_top_write_vec(t->hbond_bcoeff,p.nhbond,"HBOND_BCOEF",FLAG_F,file);
  amber_top_write_vec(t->hbond_cut,p.nhbond,"HBCUT",FLAG_F,file);
  /* graph connectivities */
  amber_top_write_vec(t->atom_types,p.natom,"AMBER_ATOM_TYPE",FLAG_S,file);
  amber_top_write_vec(t->tree_chain_class,p.natom,"TREE_CHAIN_CLASSIFICATION",
    FLAG_S,file);
  amber_top_write_vec(t->join_array,p.natom,"JOIN_ARRAY",FLAG_I,file);
  /* rotations */
  amber_top_write_vec(t->irotat,p.natom,"IROTAT",FLAG_I,file);
  /* solvent */
  if (t->pointers[AMBER_POINTER_BOX]>0) {
    amber_top_write_vec(t->box_pointer,3,"SOLVENT_POINTERS",FLAG_I,file);
    amber_top_write_vec(t->box_natoms,t->box_pointer[AMBER_BOX_NMOL],
      "ATOMS_PER_MOLECULE",FLAG_I,file);
    amber_top_write_vec(t->box_params,4,"BOX_DIMENSIONS",FLAG_F,file);
    }
  /* radius set */
  amber_top_write_val(t->radius_set,"RADIUS_SET",FLAG_S,file);
  /* capping */
  if (t->pointers[AMBER_POINTER_CAP]>0)
    msg_error("I don't know how to treat file with CAP info, sorry...",1);
  /* perturbed potential */
  if (t->pointers[AMBER_POINTER_PERT]>0) 
    msg_error("I don't know how to treat file with PERT info, sorry...",1);
  if (t->pointers[AMBER_POINTER_POL]>0) 
    msg_error("I don't know how to treat file with POL info, sorry...",1);
  /* atom radii and screening */
  amber_top_write_vec(t->radii,p.natom,"RADII",FLAG_F,file);
  amber_top_write_vec(t->screen,p.natom,"SCREEN",FLAG_F,file);
  /* polarizabilities */
  if (t->polar)
    amber_top_write_vec(t->polar,p.natom,"POLARIZABILITY",FLAG_F,file);
  amber_top_write_val(&(t->pol),"IPOL",FLAG_I,file);
  /* close the file */
  if (file_name)
    file_close(file);
  }
/* -------------------------------------------------------------------------- */
