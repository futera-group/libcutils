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
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* read force-field parameters of one atom from the file

   d  - the molecular force-field data
   id - the atom ID
   f  - open file stream */
void dlpoly_fld_read_atom(struct dlpoly_mol *d, unsigned id, FILE *f) {
  char *line,**s;
  unsigned n;
  /* atom specification */
  line =str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read atom #%d specification",1,id+1);
  s = str_split(line,' ',&n);
  if (n<3)
    msg_error_f("invalid specification of #%d atom",1,id+1); 
  line = str_free(line);
  /* name, mass, charge */
  d->atom[id].name = str_copy_new(s[0]);
  if (sscanf(s[1],"%lf",&(d->atom[id].mass))!=1)
    msg_error_f("invalid mass specification of #%d atom",1,id+1);
  if (sscanf(s[2],"%lf",&(d->atom[id].charge))!=1)
    msg_error_f("invalid charge specification of #%d atom",1,id+1);
  /* repeating counter */
  if (n>3) {
    if (sscanf(s[3],"%u",&(d->atom[id].n_atoms))!=1)
      msg_error_f("invalid counter specification of #%d atom",1,id+1);
    }
  else
    d->atom[id].n_atoms = 1;
  /* constrain */
  if (n>4) {
    if (sscanf(s[4],"%hd",&(d->atom[id].frozen))!=1)
      msg_error_f("invalid constrain specification of #%d atom",1,id+1);
    }
  else
    d->atom[id].frozen = 0;
  /* clean memory */
  vec_sfree(s,n);
  }

/* read force-field parameters of one interatomic bond from the file

   d  - the molecular force-field data
   id - the bond ID
   f  - open file stream */
void dlpoly_fld_read_bond(struct dlpoly_mol *d, unsigned id, FILE *f) {
  char *line,**s;
  unsigned n;
  /* bond specification */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read bond #%d specification",1,id+1);
  s = str_split(line,' ',&n);
  if (n<3)
    msg_error_f("invalid format of bond #%d specification",1,id+1);
  line = str_free(line);
  /* potential */
  d->bond[id].pot = dlpoly_prm_bond_pot_id(s[0]);
  /* atom IDs */
  if (sscanf(s[1],"%u",&(d->bond[id].id1))!=1)
    msg_error_f("invalid format of bond #%d ID1",1,id+1);
  if (d->bond[id].id1<1 || d->bond[id].id1>d->n_atoms)
    msg_error_f("invalid ID1 (%d) of in bond %d specification",1,
     d->bond[id].id1,id+1);
  d->bond[id].id1--;
  if (sscanf(s[2],"%u",&(d->bond[id].id2))!=1)
    msg_error_f("invalid format of bond #%d ID2",1,id+1);
  if (d->bond[id].id2<1 || d->bond[id].id2>d->n_atoms)
    msg_error_f("invalid ID2 (%d) of in bond %d specification",1,
     d->bond[id].id2,id+1);
  d->bond[id].id2--;
  /* parameters */
  if (n>3 && sscanf(s[3],"%lf",&(d->bond[id].v1))!=1)
    msg_error_f("invalid format of bond #%d parameter #1",1,id+1);
  if (n>4 && sscanf(s[4],"%lf",&(d->bond[id].v2))!=1)
    msg_error_f("invalid format of bond #%d parameter #2",1,id+1);
  if (n>5 && sscanf(s[5],"%lf",&(d->bond[id].v3))!=1)
    msg_error_f("invalid format of bond #%d parameter #3",1,id+1);
  if (n>6 && sscanf(s[6],"%lf",&(d->bond[id].v4))!=1)
    msg_error_f("invalid format of bond #%d parameter #4",1,id+1);
  /* clean memory */
  vec_sfree(s,n);
  }

/* read force-field parameters of one valence angle from the file

   d  - the molecular force-field data
   id - the angle ID
   f  - open file stream */
void dlpoly_fld_read_angle(struct dlpoly_mol *d, unsigned id, FILE *f) {
  char *line,**s;
  unsigned n;
  /* angle specification */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read angle #%d specification",1,id+1);
  s = str_split(line,' ',&n);
  if (n<4)
    msg_error_f("invalid format of angle #%d specification",1,id+1);
  line = str_free(line);
  /* potential */
  d->angle[id].pot = dlpoly_prm_angle_pot_id(s[0]);
  /* atom IDs */
  if (sscanf(s[1],"%u",&(d->angle[id].id1))!=1)
    msg_error_f("invalid format of angle #%d ID1",1,id+1);
  if (d->angle[id].id1<1 || d->angle[id].id1>d->n_atoms)
    msg_error_f("invalid ID1 (%d) of in angle %d specification",1,
     d->angle[id].id1,id+1);
  d->angle[id].id1--;
  if (sscanf(s[2],"%u",&(d->angle[id].id2))!=1)
    msg_error_f("invalid format of angle #%d ID2",1,id+1);
  if (d->angle[id].id2<1 || d->angle[id].id2>d->n_atoms)
    msg_error_f("invalid ID2 (%d) of in angle %d specification",1,
     d->angle[id].id2,id+1);
  d->angle[id].id2--;
  if (sscanf(s[3],"%u",&(d->angle[id].id3))!=1)
    msg_error_f("invalid format of angle #%d ID3",1,id+1);
  if (d->angle[id].id3<1 || d->angle[id].id3>d->n_atoms)
    msg_error_f("invalid ID3 (%d) of in angle %d specification",1,
     d->angle[id].id3,id+1);
  d->angle[id].id3--;
  /* parameters */
  if (n>4 && sscanf(s[4],"%lf",&(d->angle[id].v1))!=1)
    msg_error_f("invalid format of angle #%d parameter #1",1,id+1);
  if (n>5 && sscanf(s[5],"%lf",&(d->angle[id].v2))!=1)
    msg_error_f("invalid format of angle #%d parameter #2",1,id+1);
  if (n>6 && sscanf(s[6],"%lf",&(d->angle[id].v3))!=1)
    msg_error_f("invalid format of angle #%d parameter #3",1,id+1);
  if (n>7 && sscanf(s[7],"%lf",&(d->angle[id].v4))!=1)
    msg_error_f("invalid format of angle #%d parameter #4",1,id+1);
  /* clean memory */
  vec_sfree(s,n);
  }

/* read force-field parameters of one dihedral angle from the file

   d  - the molecular force-field data
   id - the dihedral ID
   f  - open file stream */
void dlpoly_fld_read_dihed(struct dlpoly_mol *d, unsigned id, FILE *f) {
  char *line,**s;
  unsigned n;
  /* dihedral angle specification */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read dihedral angle #%d specification",1,id+1);
  s = str_split(line,' ',&n);
  if (n<5)
    msg_error_f("invalid format of dihedral angle #%d specification",1,id+1);
  line = str_free(line);
  /* potential */
  d->dihed[id].pot = dlpoly_prm_dihed_pot_id(s[0]);
  /* atom IDs */
  if (sscanf(s[1],"%u",&(d->dihed[id].id1))!=1)
    msg_error_f("invalid format of dihedral angle #%d ID1",1,id+1);
  if (d->dihed[id].id1<1 || d->dihed[id].id1>d->n_atoms)
    msg_error_f("invalid ID1 (%d) of in dihedral angle %d specification",1,
     d->dihed[id].id1,id+1);
  d->dihed[id].id1--;
  if (sscanf(s[2],"%u",&(d->dihed[id].id2))!=1)
    msg_error_f("invalid format of dihedral angle #%d ID2",1,id+1);
  if (d->dihed[id].id2<1 || d->dihed[id].id2>d->n_atoms)
    msg_error_f("invalid ID2 (%d) of in dihedral angle %d specification",1,
     d->dihed[id].id2,id+1);
  d->dihed[id].id2--;
  if (sscanf(s[3],"%u",&(d->dihed[id].id3))!=1)
    msg_error_f("invalid format of dihedral angle #%d ID3",1,id+1);
  if (d->dihed[id].id3<1 || d->dihed[id].id3>d->n_atoms)
    msg_error_f("invalid ID3 (%d) of in dihedral angle %d specification",1,
     d->dihed[id].id3,id+1);
  d->dihed[id].id3--;
  if (sscanf(s[4],"%u",&(d->dihed[id].id4))!=1)
    msg_error_f("invalid format of dihedral angle #%d ID4",1,id+1);
  if (d->dihed[id].id4<1 || d->dihed[id].id4>d->n_atoms)
    msg_error_f("invalid ID4 (%d) of in dihedral angle %d specification",1,
     d->dihed[id].id4,id+1);
  d->dihed[id].id4--;
  /* parameters */
  if (n> 5 && sscanf(s[ 5],"%lf",&(d->dihed[id].v1))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #1",1,id+1);
  if (n> 6 && sscanf(s[ 6],"%lf",&(d->dihed[id].v2))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #2",1,id+1);
  if (n> 7 && sscanf(s[ 7],"%lf",&(d->dihed[id].v3))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #3",1,id+1);
  if (n> 8 && sscanf(s[ 8],"%lf",&(d->dihed[id].v4))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #4",1,id+1);
  if (n> 9 && sscanf(s[ 9],"%lf",&(d->dihed[id].v5))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #5",1,id+1);
  if (n>10 && sscanf(s[10],"%lf",&(d->dihed[id].v6))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #6",1,id+1);
  if (n>11 && sscanf(s[11],"%lf",&(d->dihed[id].v7))!=1)
    msg_error_f("invalid format of dihedral angle #%d parameter #7",1,id+1);
  /* clean memory */
  vec_sfree(s,n);
  }

/* read force-field parameters of one molecule from the file

   d  - the force-field data
   id - the molecule ID
   f  - open file stream */
void dlpoly_fld_read_mol(struct dlpoly_fld *d, unsigned id, FILE *f) {
  char *line,key[1024];
  unsigned i,n;
  /* name of the molecule */
  d->mol[id].name = str_read_line_new(f);
  if (!d->mol[id].name)
    msg_error_f("cannot read name of molecule #%d",1,id+1);
  str_trim(d->mol[id].name);
  /* number of molecules */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read number of #%d molecules",1,id+1);
  if (sscanf(line,"%s%u",key,&(d->mol[id].n_molecules))!=2)
    msg_error_f("invalid number of #%d molecules specification",1,id+1);
  str_lowcase(key);
  if (!str_compare("nummols",key))
    msg_error_f("NUMMOLS keyword expected (\"%s\" read)",1,key);
  /* number of atoms */
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read number of #%d atoms",1,id+1);
  if (sscanf(line,"%s%u",key,&(d->mol[id].n_atoms))!=2)
    msg_error_f("invalid number of #%d atoms specification",1,id+1);
  str_lowcase(key);
  if (!str_compare("atoms",key))
    msg_error_f("ATOMS keyword expected (\"%s\" read)",1,key);
  line = str_free(line);
  d->mol[id].atom = dlpoly_atom_new(d->mol[id].n_atoms);
  for (i=0, n=0; i<d->mol[id].n_atoms && n<d->mol[id].n_atoms; i++) {
    dlpoly_fld_read_atom(&(d->mol[id]),i,f);
    n = n + d->mol[id].atom[i].n_atoms;
    }
  /* force-field parameters */
  for (line=str_read_line_new(f); line;
       line=str_free(line), line=str_read_line_new(f)) {
    /* keyword */
    if (sscanf(line,"%s",key)!=1)
      msg_error_f("molecular keyword expected (\"%s\" read)",1,key);
    str_lowcase(key);
    /* core shells */
    if (str_compare(key,"shell"))
      msg_error("shell potentials are not supported",1);
    /* constraints */
    else if (str_compare(key,"constraints"))
      msg_error("constrain potentials are not supported",1);
    /* PMF bond lengths */
    else if (str_compare(key,"pmf"))
      msg_error("pmf potentials are not supported",1);
    /* rigid units */
    else if (str_compare(key,"rigid"))
      msg_error("rigid units are not supported",1);
    /* tethered atoms */
    else if (str_compare(key,"teth"))
      msg_error("tethering potentials are not supported",1);
    /* interatomic bonds */
    else if (str_compare(key,"bonds")) {
      if (sscanf(line,"%s%u",key,&(d->mol[id].n_bonds))!=2)
        msg_error("invalid format of number-of-bonds specification",1);
      d->mol[id].bond = dlpoly_prm_bond_new(d->mol[id].n_bonds);
      for (i=0; i<d->mol[id].n_bonds; i++) 
        dlpoly_fld_read_bond(&(d->mol[id]),i,f);
      }
    /* valence angles */
    else if (str_compare(key,"angles")) {
      if (sscanf(line,"%s%u",key,&(d->mol[id].n_angles))!=2)
        msg_error("invalid format of number-of-angles specification",1);
      d->mol[id].angle = dlpoly_prm_angle_new(d->mol[id].n_angles);
      for (i=0; i<d->mol[id].n_angles; i++) 
        dlpoly_fld_read_angle(&(d->mol[id]),i,f);
      }
    /* dihedral angles */
    else if (str_compare(key,"dihedrals")) {
      if (sscanf(line,"%s%u",key,&(d->mol[id].n_dihedrals))!=2)
        msg_error("invalid format of number-of-dihedrals specification",1);
      d->mol[id].dihed = dlpoly_prm_dihed_new(d->mol[id].n_dihedrals);
      for (i=0; i<d->mol[id].n_dihedrals; i++) 
        dlpoly_fld_read_dihed(&(d->mol[id]),i,f);
      }
    /* inversions */
    else if (str_compare(key,"inversions"))
      msg_error("inversion potentials are not supported",1);
    /* end of the molecular specification */
    else if (str_compare(key,"finish"))
      break;
    else
      msg_error_f("invalid molecular keyword \"%s\"",1,key);
    }
  }

/* read force-field parameters of one VDW potential from the file

   d  - the force-field data
   id - the potential ID
   f  - open file stream */
void dlpoly_fld_read_vdw(struct dlpoly_fld *d, unsigned id, FILE *f) {
  char *line,**s;
  unsigned i,j,n;
  /* potential specification */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read VDW potential #%d specification",1,id+1);
  s = str_split(line,' ',&n);
  if (n<3)
    msg_error_f("invalid format of VDW potential #%d specification",1,id+1);
  line = str_free(line);
  /* potential */
  d->vdw[id].pot = dlpoly_prm_vdw_pot_id(s[2]);
  /* atom types */
  d->vdw[id].at1 = str_copy_new(s[0]);
  if (!dlpoly_fld_atom_find(d,d->vdw[id].at1,&i,&j))
    msg_error_f("unknown atom type \"%s\" in VDW potential #%d",
      1,d->vdw[id].at1,id+1);
  d->vdw[id].at2 = str_copy_new(s[1]);
  if (!dlpoly_fld_atom_find(d,d->vdw[id].at2,&i,&j))
    msg_error_f("unknown atom type \"%s\" in VDW potential #%d",
      1,d->vdw[id].at2,id+1);
  /* parameters */
  if (n>3 && sscanf(s[3],"%lf",&(d->vdw[id].v1))!=1)
    msg_error_f("invalid format of VDW potential #%d parameter #1",1,id+1);
  if (n>4 && sscanf(s[4],"%lf",&(d->vdw[id].v2))!=1)
    msg_error_f("invalid format of VDW potential #%d parameter #2",1,id+1);
  if (n>5 && sscanf(s[5],"%lf",&(d->vdw[id].v3))!=1)
    msg_error_f("invalid format of VDW potential #%d parameter #3",1,id+1);
  if (n>6 && sscanf(s[6],"%lf",&(d->vdw[id].v4))!=1)
    msg_error_f("invalid format of VDW potential #%d parameter #4",1,id+1);
  if (n>7 && sscanf(s[7],"%lf",&(d->vdw[id].v5))!=1)
    msg_error_f("invalid format of VDW potential #%d parameter #5",1,id+1);
  /* clean memory */
  vec_sfree(s,n);
  }

/* read force-field parameters from the external file
 
   d    - the force-field data
   name - name of the file */
void dlpoly_fld_read(struct dlpoly_fld *d, char *name) {
  char *line,key[1024],val[1024];
  unsigned i;
  FILE *f;
  f = file_open(name,"r");
  /* header */
  d->header = str_read_line_new(f);
  if (!d->header)
    msg_error_f("cannot read header from \"%s\" file",1,name);
  str_trim(d->header);
  /* units */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read units from \"%s\" file",1,name);
  if (sscanf(line,"%s%s",key,val)!=2)
    msg_error("invalid format of unit specification",1);
  str_lowcase(key);
  if (!str_compare("units",key))
    msg_error_f("UNITS keyword expected (\"%s\" read)",1,key);
  d->units = dlpoly_fld_unit_id(val);
  line = str_free(line);
  /* force-field data */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    /* ignore blank lines */
    if (str_length(line)<2)
      continue;
    /* keyword */
    if (sscanf(line,"%s",key)!=1)
      msg_error_f("system keyword expected (\"%s\" read)",1,key);
    str_lowcase(key);
    /* molecular data */
    if (str_compare(key,"molecules")) {
      if (sscanf(line,"%s%u",key,&(d->n_molecules))!=2)
        msg_error("invalid format of number-of-molecule specification",1);
      d->mol = dlpoly_mol_new(d->n_molecules);
      for (i=0; i<d->n_molecules; i++) 
        dlpoly_fld_read_mol(d,i,f);
      }
    /* OpenKIM model */
    else if (str_compare(key,"kim"))
      msg_error("OpenKIM model is not supported",1);
    /* VdW pair potentials */
    else if (str_compare(key,"vdw")) {
      if (sscanf(line,"%s%u",key,&(d->n_vdws))!=2)
        msg_error("invalid format of number-of-vdw-pot specification",1);
      d->vdw = dlpoly_prm_vdw_new(d->n_vdws);
      for (i=0; i<d->n_vdws; i++)
        dlpoly_fld_read_vdw(d,i,f);
      }
    /* metal potentials */
    else if (str_compare(key,"metal"))
      msg_error("metal potentials are not supported",1);
    /* RDF pairs */
    else if (str_compare(key,"rdf"))
      msg_error("RDF pairs are not supported",1);
    /* tersoff potentials */
    else if (str_compare(key,"tersoff"))
      msg_error("tersoff potentials are not supported",1);
    /* three-body potentials */
    else if (str_compare(key,"tbp"))
      msg_error("three-body potentials are not supported",1);
    /* four-body potentials */
    else if (str_compare(key,"fbp"))
      msg_error("four-body potentials are not supported",1);
    /* external fields */
    else if (str_compare(key,"extern"))
      msg_error("external field are not supported",1);
    /* end of the file */
    else if (str_compare(key,"close"))
      break;
    /* unknown keyword */
    else
      msg_error_f("invalid keyword \"%s\" found in \"%s\" file",1,key,name);
    /* end of the file */
    }
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
