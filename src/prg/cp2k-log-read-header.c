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

#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read header data from cp2k output file
 
   d - pointer to cp2k data struct
   f - open file stream */
void cp2k_log_read_header(struct cp2k_dat *d, FILE *f) {
  char *line;
  unsigned i,n;
  double cv[3][3];
  fpos_t fpos;
  /* version */
  line = str_ffind_b_new(f," CP2K| version string:");
  if (!line)
    msg_error("cp2k version specification not found",1);
  str_trim(line);
  for (n=str_length(line); n>0 && line[n-1]!=' '; n--);
  d->version = str_copy_new(line+n);
  str_free(line);
  /* project name */
  line = str_ffind_b_new(f," GLOBAL| Project name");
  if (!line)
    msg_error("cp2k global project name not found",1);
  str_trim(line);
  for (n=str_length(line); n>0 && line[n-1]!=' '; n--);
  d->project_name = str_copy_new(line+n);
  str_free(line);
  /* run type */
  line = str_ffind_b_new(f," GLOBAL| Run type");
  if (!line)
    msg_error("cp2k global run type not found",1);
  str_trim(line);
  for (n=str_length(line); n>0 && line[n-1]!=' '; n--);
  d->run_type = cp2k_run_type_id(line+n);
  str_free(line);
  /* print level */
  line = str_ffind_b_new(f," GLOBAL| Global print level");
  if (!line)
    msg_error("cp2k global print level not found",1);
  str_trim(line);
  for (n=str_length(line); n>0 && line[n-1]!=' '; n--);
  d->print_level = cp2k_print_level_id(line+n);
  str_free(line);
  /* cell */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," CELL| Volume");
  if (line) {
    line = str_free(line);
    for (i=0; i<3; i++) {
      line = str_read_line_new(f); 
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      if (sscanf(line+28,"%lf%lf%lf",&(cv[i][0]),&(cv[i][1]),&(cv[i][2]))!=3)
        msg_error_f("invalid format of cell vector #%d",1,i+1);
      line = str_free(line);
      }
    cell_set_vec(d->cell,cv[0],cv[1],cv[2]);
    }
  else
    fsetpos(f,&fpos);
  /* K-points */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," BRILLOUIN| List of Kpoints");
  if (line) {
    if (sscanf(line+50,"%d",&(d->n_kpoints))!=1)
      msg_error("cannot read number of k-points",1);
    line = str_free(line);
    line = str_read_line_new(f);
    line = str_free(line);
    d->kpoint = cp2k_kpoint_new(d->n_kpoints);
    for (i=0; i<d->n_kpoints; i++) {
      line = str_read_line_new(f);
      if (sscanf(line+20,"%lf%lf%lf%lf",&(d->kpoint[i].weight),
           &(d->kpoint[i].crd[0]),&(d->kpoint[i].crd[1]),
           &(d->kpoint[i].crd[2]))!=4)
        msg_error_f("cannot read k-point #%d coordinates and weight",1,i+1);
      line = str_free(line);
      }
    }
  else
    fsetpos(f,&fpos);
  /* spin multiplicity */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," DFT| Multiplicity");
  if (line) {
    if (sscanf(line+30,"%u",&(d->multiplicity))!=1)
      msg_error("cannot read spin multiplicity",1);
    line = str_free(line);
    }
  else {
    d->multiplicity = 1;
    fsetpos(f,&fpos);
    }
  /* number of spin components */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," DFT| Number of spin states");
  if (line) {
    if (sscanf(line+30,"%u",&(d->n_spins))!=1)
      msg_error("cannot read number of spin components",1);
    line = str_free(line);
    }
  else {
    fsetpos(f,&fpos);
    line = str_ffind_b_new(f," Spin 2");
    if (line)
      d->n_spins = 2;
    else
      d->n_spins = 1;
    fsetpos(f,&fpos);
    }
  /* net charge */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," DFT| Charge");
  if (line) {
    if (sscanf(line+30,"%u",&(d->charge))!=1)
      msg_error("cannot read net charge",1);
    line = str_free(line);
    }
  else {
    d->charge = 0;
    fsetpos(f,&fpos);
    }
  /* number of atoms and orbitals */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," TOTAL NUMBERS AND MAXIMUM NUMBERS");
  if (line) {
    line = str_free(line);
    for (i=0; i<8; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      if (i>0) {
        if (sscanf(line+70,"%d",&n)!=1)
          msg_error("cannot read data in total/maximum number block",1);
        switch (i) {
          case 1: d->n_atom_kinds = n; break;
          case 2: d->n_atoms = n;      break;
          case 3: d->n_shell_sets = n; break;
          case 4: d->n_shells = n;     break;
          case 5: d->n_prim_fce = n;   break;
          case 6: d->n_cart_fce = n;   break;
          case 7: d->n_sphr_fce = n;   break;
          }
        line = str_free(line);
        }
      }
    }
  else
    fsetpos(f,&fpos);
  /* number of electrons / states */
  d->n_electrons = vec_ualloc(d->n_spins);
  vec_uset(d->n_electrons,0,d->n_spins);
  d->n_states = vec_ualloc(d->n_spins);
  vec_uset(d->n_states,0,d->n_spins);
  if (d->n_spins>1) {
    /* alpha spin */
    fgetpos(f,&fpos);
    line = str_ffind_b_new(f," Spin 1");  
    if (line) {
      line = str_free(line);
      line = str_ffind_b_new(f," Number of electrons:");
      if (sscanf(line+70,"%d",&(d->n_electrons[0]))!=1)
        msg_error("cannot read number of alpha electrons",1);
      line = str_free(line);
      line = str_ffind_b_new(f," Number of molecular orbitals:");
      if (sscanf(line+70,"%d",&(d->n_states[0]))!=1)
        msg_error("cannot read number of alpha states",1);
      line = str_free(line);
      }
    else
      fsetpos(f,&fpos);
    /* beta spin */
    fgetpos(f,&fpos);
    line = str_ffind_b_new(f," Spin 2");  
    if (line) {
      line = str_free(line);
      line = str_ffind_b_new(f," Number of electrons:");
      if (sscanf(line+70,"%d",&(d->n_electrons[1]))!=1)
        msg_error("cannot read number of beta electrons",1);
      line = str_free(line);
      line = str_ffind_b_new(f," Number of molecular orbitals:");
      if (sscanf(line+70,"%d",&(d->n_states[1]))!=1)
        msg_error("cannot read number of beta states",1);
      line = str_free(line);
      }
    else
      fsetpos(f,&fpos);
    }
  else {
    /* electrons */
    fgetpos(f,&fpos);
    line = str_ffind_b_new(f," Number of electrons:");
    if (line) {
      if (sscanf(line+70,"%d",&(d->n_electrons[0]))!=1)
        msg_error("cannot read number of electrons",1);
      }
    else
      fsetpos(f,&fpos);
    line = str_free(line);
    /* states */
    fgetpos(f,&fpos);
    line = str_ffind_b_new(f," Number of molecular orbitals:");
    if (line) {
      if (sscanf(line+70,"%d",&(d->n_states[0]))!=1)
        msg_error("cannot read number of states",1);
      }
    else
      fsetpos(f,&fpos);
    line = str_free(line);
    }
  /* number of basis functions (AOs) */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," Number of orbital functions:");
  if (line) {
    if (sscanf(line+70,"%d",&(d->n_bfce))!=1)
      msg_error("cannot read number of orbital functions",1);
    }
  else
    fsetpos(f,&fpos);
  line = str_free(line);
  /* number of independent basis functions (AOs) */
  fgetpos(f,&fpos);
  line = str_ffind_b_new(f," Number of independent orbital functions:");
  if (line) {
    if (sscanf(line+70,"%d",&(d->n_ibfce))!=1)
      msg_error("cannot read number of independent orbital functions",1);
    }
  else
    fsetpos(f,&fpos);
  line = str_free(line);
  if (!d->n_ibfce)
    d->n_ibfce = d->n_sphr_fce;
  }

/* -------------------------------------------------------------------------- */
