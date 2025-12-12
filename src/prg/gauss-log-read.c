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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <cmn/file.h>
#include <cmn/list.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include <cmn/units.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* read one line from gaussian log file
 
   w - descpription of what is reading
   f - pointer to open gaussian log file */
char* gauss_log_read_line(char *w, FILE *f) {
  char *line;
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("unexpected end of gaussian log file while reading %s",1,w);
  return(line);
  }

/* return internal code of computational method

   g - pointer to gaussian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_method(struct gauss_dat *g, FILE *f) {
  char *line;
  unsigned n;
  /* list of keywords */
  line = gauss_log_read_line("method specification",f);
  n = str_length(line);
  if (n<5 || line[1]!='#')
    msg_error("invalid format of control line",1);
  g->job_mthd_id = gauss_mthd_id(line);
  line = str_free(line);
  /* list of IOps */
  line = gauss_log_read_line("method specification footline",f);
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    n = str_length(line);
    if (n<3 || line[n-2]!=';') {
      line = str_free(line);
      break;
      }
    gauss_iop_set(g->iop,line);
    }
  }

/* read main info (program code, method, basis set)

   g - pointer to gaussian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_main(struct gauss_dat *g, FILE *f) {
  char program[80],revision[80];
  char *line = NULL;
  short prg_found = 0,mthd_found = 0;
  rewind(f);
  /* read molecule specification */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    /* program specification */
    if (!prg_found &&
        strstr(line," ******************************************")) {
      line = gauss_log_read_line("program code specification",f);
      if (sscanf(line,"%*s%s%s",program,revision)!=2)
        msg_error("cannot recognize gaussian file format",1);
      g->prog_ver = gauss_prg_ver_id(program);
      g->prog_rev = gauss_prg_rev_id(revision);
      prg_found = 1;
      line = str_free(line);
      }
    /* method specification */
    if (prg_found && !mthd_found &&
        line && strstr(line," ---")) {
      gauss_log_read_method(g,f);
      mthd_found = 1;
      }
    if (prg_found && mthd_found)
      break;
    }
  }

/* read molecule specification from open gaussian log file

   g - pointer to gaussian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_spec_f(struct gauss_dat *g, FILE *f) {
  char **w,*line = NULL;
  unsigned n;
  short found = 0;
  /* read main info */
  gauss_log_read_main(g,f);
  /* read molecule specification */
  rewind(f);
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f))
    if (strstr(line,"Charge =")) {
      w = str_split(line,' ',&n);
      if (n!=6)
        msg_error("invalid format of charge and multipliticy",1);
      type_read(w[2],&(g->charge),TYPE_INT);
      type_read(w[5],&(g->multiplicity),TYPE_UINT);
      vec_sfree(w,n);
      }
    else if (strstr(line,"NAtoms=")) {
      if (sscanf(line+8,"%d",&(g->n_atoms))!=1)
        msg_error("cannot read number of atom specification in input file",1);
      found = 1;
      str_free(line);
      break;
      }
  /* check if required dat was found */
  if (!found)
    msg_error("cannot find molecule specification data in input file",1);
  found = 0;
  /* read number of degrees of freedom */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f))
    if (strstr(line,"Deg. of freedom")) {
      if (sscanf(line+16,"%d",&(g->n_fdeg))!=1)
        msg_error("cannot read number of degrees of freedom in input file",1);
      found = 1;
      str_free(line);
      break;
      }
  /* check if required data was found */
  if (!found)
    msg_error("cannot find number of degrees of freedom in input file",1);
  /* number of vibrational modes */
  g->n_vib_modes = (3*g->n_atoms)-6;
  }

/* read molecule specification from gaussian log file

   g - pointer to gaussian data struct
   f - name of the log file */
void gauss_log_read_spec(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_spec_f(g,file);
  file_close(file);
  }

/* read NMR shifts from open gaussian log file
 
   g - pointer to guassian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_nmr_f(struct gauss_dat *g, FILE *f) {
  char *line,msg[80];
  unsigned ia,ic;
  rewind(f);
  sprintf(msg,"unexpected end of file while reading nmr shifts");
  /* read nmr shifts */
  line = str_ffind_new(f,"SCF GIAO Magnetic shielding tensor (ppm):");
  if (line) {
    line = str_free(line);
    /* reallocate memory */
    g->nmr = gauss_nmr_free(g->nmr,g->n_atoms);
    g->nmr = gauss_nmr_new(g->n_atoms);
    /* read magnetic shield tensors */
    for (ia=0; ia<g->n_atoms; ia++) {
      /* isotropic and anisotropic part */
      line = str_read_line_new(f);
      if (!line)
        msg_error(msg,1);
      if (sscanf(line,"%*d%*s%*s%*s%lf%*s%*s%lf",
          &(g->nmr[ia].isotropic),&(g->nmr[ia].anisotropy))!=2)
        msg_error("cannot read mag. shield tensor elements",1);
      line = str_free(line);
      /* mag. shield tensor */
      for (ic=0; ic<3; ic++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error(msg,1);
        if (sscanf(line+7,"%lf",&(g->nmr[ia].tensor[ic][0]))!=1)
          msg_error("cannot read x-component of mag. shield tensor",1);
        if (sscanf(line+24,"%lf",&(g->nmr[ia].tensor[ic][1]))!=1)
          msg_error("cannot read y-component of mag. shield tensor",1);
        if (sscanf(line+41,"%lf",&(g->nmr[ia].tensor[ic][2]))!=1)
          msg_error("cannot read z-component of mag. shield tensor",1);
        line = str_free(line);
        }
      /* eigenvalues */
      line = str_read_line_new(f);
      if (!line)
        msg_error(msg,1);
      if (sscanf(line+16,"%lf%lf%lf",&(g->nmr[ia].eigen[0]),
        &(g->nmr[ia].eigen[1]),&(g->nmr[ia].eigen[2]))!=3)
        msg_error("cannot read eigenvalues of mag. shield tensor",1);
      }
    }
  }

/* read NMR shifts from gaussian log file
 
   g - pointer to guassian data struct
   f - name of the gaussian log file */
void gauss_log_read_nmr(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_nmr_f(g,file);
  file_close(file);
  }

/* read rotatory strength for all excited states
 
   g - pointer to guassian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_states_rot_f(struct gauss_dat *g, FILE *f) {
  short found_vel = 0,found_len = 0,process_line = 0;
  char *line,msg[80];
  unsigned i,j;
  rewind(f);
  sprintf(msg,"unexpected end of file while reading rotatory strengths");
  /* read the file */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    /* rotatory strengths */
    if (strstr(line,"Rotatory Strengths (R) in cgs")) {
      line = str_free(line);
      line = str_read_line_new(f);
      if (!line)
        msg_error(msg,1);
      /* velocity definition */
      if (strstr(line+50,"R(velocity)")) {
        found_vel = 1;
        for (i=0; i<g->n_states; i++) {
          /* strength */
          if (process_line)
            process_line = 0;
          else {
            line = str_free(line);
            line = str_read_line_new(f);
            if (!line)
              msg_error(msg,1);
            }
          if (sscanf(line+50,"%lf",&(g->state[i].rot_vel))!=1)
            msg_error_f("cannot read rot.s. (velocity) of state #%d",1,i+1);
          /* check format */
          line = str_free(line);
          line = str_read_line_new(f);
          if (!line)
            msg_error(msg,1);
          if (strstr(line,"tensor")) {
            /* total R(velocity) tensor */
            line = str_free(line);
            for (j=0; j<4; j++) {
              line = str_read_line_new(f);
              if (!line)
                msg_error(msg,1);
              line = str_free(line);
              }
            /* tensor in input orientation */
            line = str_free(line);
            for (j=0; j<5; j++) {
              line = str_read_line_new(f);
              if (!line)
                msg_error(msg,1);
              line = str_free(line);
              }
            }
          else
            process_line = 1;
          }
        }
      /* dipole length definition */
      else if (strstr(line+50,"R(length)")) {
        for (i=0; i<g->n_states; i++) {
          line = str_free(line);
          line = str_read_line_new(f);
          if (!line)
            msg_error(msg,1);
          if (sscanf(line+50,"%lf",&(g->state[i].rot_len))!=1)
            msg_error_f("cannot read rot.s. (length) of state #%d",1,i+1);
          }
        found_len = 1;
        }
      /* invalid format */
      else
        msg_error("invalid format of rotatory strengt format",1);
      if (found_vel && found_len)
        break;
      }
    }
  }

/* read excited states from open gaussian log file
 
   g - pointer to guassian data struct
   f - pointer to open gaussian log file */
void gauss_log_read_states_f(struct gauss_dat *g, FILE *f) {
  char *line,state[1024];
  unsigned i,id = 0;
  double *d,d1,d2,d3,ds,os;
  struct gauss_state *s;
  struct queue *qe = NULL,*qd = NULL,*qx = NULL;
  rewind(f);
  /* read the file */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    /* skip lower methods */
    if (g->job_mthd_id==GAUSS_MTHD_EOMCCSD) {
      line = str_free(line);
      line = str_ffind_new(f,"EOM-CCSD transition properties");
      if (!line)
        msg_warn("EOM-CCSD data reading is not supported, skipped");
      line = str_free(line);
      break;
      }
    /* read ground-to-excited-state transition dipoles */
    if (line && strstr(line,"Ground to excited state transition electric dipole moments")) {
      qd = queue_alloc();
      line = str_free(line);
      line = str_read_line_new(f);
      if (!line)
        break;
      for (line=str_free(line),line=str_read_line_new(f); line;
           line=str_free(line),line=str_read_line_new(f)) {
        if (sscanf(line,"%*d%lf%lf%lf%lf%lf",&d1,&d2,&d3,&ds,&os)!=5) {
          line = str_free(line);
          break;
          }
        d = vec_falloc(3);
        d[0] = d1;
        d[1] = d2;
        d[2] = d3;
        queue_add(qd,d);
        }
      }
    /* read excited-to-excited-state transition dipoles */
    if (line && strstr(line,"Excited to excited state transition electric dipole moments")) {
      qx = queue_alloc();
      line = str_free(line);
      line = str_read_line_new(f);
      if (!line)
        break;
      for (line=str_free(line),line=str_read_line_new(f); line;
           line=str_free(line),line=str_read_line_new(f)) {
        if (strstr(line,"States")) {
          if (sscanf(line,"%*s%*d%*s%*d%lf%lf%lf%lf",&d1,&d2,&d3,&os)!=4) {
            line = str_free(line);
            break;
            }
          }
        else if (sscanf(line,"%*d%*d%lf%lf%lf%lf%lf",&d1,&d2,&d3,&ds,&os)!=5) {
          line = str_free(line);
          break;
          }
        d = vec_falloc(3);
        d[0] = d1;
        d[1] = d2;
        d[2] = d3;
        queue_add(qx,d);
        }
      }
    /* read excitation energies to queue */
    if (line && strstr(line,"Excitation energies and oscillator strengths")) {
      qe = queue_alloc();
      for (line=str_free(line),line=str_read_line_new(f); line;
           line=str_free(line),line=str_read_line_new(f)) {
        if (strstr(line," ************")) {
          line = str_free(line);
          break;
          }
        if (strstr(line," Excited State ")) {
          s = gauss_state_new(1);
          if (sscanf(line+20,"%s%lf eV %lf nm  f=%lf <S**2>=%lf",
            state,&(s->energy),&(s->lambda),&(s->strength),&(s->spin_s2))!=5)
            msg_error("cannot read excited state data (general format)",1);
          queue_add(qe,s);
          }
        }
      break;
      }
    }
  /* convert queue to data array */
  if (qe) {
    if (qd && qe->num!=qd->num)
      msg_error_f("inconsistent number of excited state (%ld)"
        " and transition dipoles (%ld)",1,qe->num,qd->num);
    if (qx && qx->num!=(qe->num*qe->num-qe->num)/2)
      msg_error_f("%ld excited-to-excited-state transition dipoles"
        " found, %ld expected",1,qx->num,(qe->num*qe->num-qe->num)/2);
    g->n_states = qe->num;
    g->state = gauss_state_new(g->n_states);
    while (qe->num) {
      s = queue_get(qe);
      g->state[id].energy = s->energy;
      g->state[id].lambda = s->lambda;
      g->state[id].strength = s->strength;
      g->state[id].spin_s2 = s->spin_s2;
      gauss_state_free(s,1);
      if (qd) {
        d = queue_get(qd);
        g->state[id].dipole = vec_fcopy_new(d,3);
        d = vec_ffree(d);
        }
      if (qx && id) {
        g->state[id].dipole_ee = vec_talloc(sizeof(double*),id);
        for (i=0; i<id; i++) {
          d = queue_get(qx);
          g->state[id].dipole_ee[i] = vec_fcopy_new(d,3);
          d = vec_ffree(d);
          }
        }
      id++;
      }
    }
  /* clean memory */
  queue_free(qe);
  queue_free(qd);
  queue_free(qx);
  /* post-HF corrections */
  id = 0;
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (strstr(line,"CIS(D) Correction to the")) {
      line = str_free(line);
      line = str_ffind_new(f,"CIS(D) Exc. E:");
      if (!line)
        msg_error_f("cannot find CIS(D) correction to state %d",1,id+1);
      if (sscanf(line+40,"%lf eV %lf nm",
          &(g->state[id].energy),&(g->state[id].lambda))!=2)
        msg_error_f("cannot read CIS(D) correction to state %d",1,id+1);
      id++;
      }
    if (id==g->n_states)
      break;
    }
   /* rotatory strengths */
  if (g->n_states)
    gauss_log_read_states_rot_f(g,f);
  }

/* read excited states from gaussian log file
 
   g - pointer to guassian data struct
   f - name of the gaussian log file */
void gauss_log_read_states(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_states_f(g,file);
  file_close(file);
  }

/* read coordinates from open gaussian log file
 
   g - guassian data struct
   c - specify coordinate orientation
   f - open gaussian-log file stream */
void gauss_log_read_coord_f(struct gauss_dat *g, short c, FILE *f) {
  char *line,*chp;
  unsigned i,j,center,type;
  double **crd,energy;
  struct queue *qs,*qe;
  qs = queue_alloc();
  qe = queue_alloc();
  /* read coordinates */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (((c==GAUSS_CRD_INP || c==GAUSS_CRD_ANY)
          && strstr(line,"Input orientation:")) ||
        ((c==GAUSS_CRD_STD || c==GAUSS_CRD_ANY)
          && strstr(line,"Standard orientation:")) ||
        ((c==GAUSS_CRD_ZMT || c==GAUSS_CRD_ANY)
          && strstr(line,"Z-Matrix orientation:"))) {
      /* skip header */
      for (i=0; i<4; i++) {
        line = str_free(line);
        line = gauss_log_read_line("coordinates",f);
        }
      /* allocate memory */
      crd = mat_falloc(g->n_atoms,3);
      if (!g->atom)
        g->atom = gauss_atom_new(g->n_atoms);
      /* coordinates */
      for (i=0; i<g->n_atoms; i++) {
        line = str_free(line);
        line = gauss_log_read_line("coordinates",f);
        if (sscanf(line,"%d%d%d%lf%lf%lf",&center,&(g->atom[i].num),&type,
          &(crd[i][0]),&(crd[i][1]),&(crd[i][2]))!=6)
          msg_error("invalid coordinate format in input file",1);
        }
      queue_add(qs,crd);
      }
    if ((c==GAUSS_CRD_TRJ || c==GAUSS_CRD_ANY)
          && strstr(line,"  Cartesian coordinates:")) {
      /* allocate memory */
      crd = mat_falloc(g->n_atoms,3);
      if (!g->atom)
        g->atom = gauss_atom_new(g->n_atoms);
      /* coordinates */
      for (i=0; i<g->n_atoms; i++) {
        line = str_free(line);
        line = gauss_log_read_line("coordinates",f);
        for (j=0; j<str_length(line); j++)
          if (line[j]=='D')
            line[j] = 'E';
        if (sscanf(line,"%*s%*s%*s%lf%*s%lf%*s%lf",
          &(crd[i][0]),&(crd[i][1]),&(crd[i][2]))!=3)
          msg_error("invalid coordinate format in input file",1);
        }
      queue_add(qs,crd);
      }
    /* energy */
    if (strstr(line,"SCF Done:")) {
      chp = strstr(line,"=");
      sscanf(chp+1,"%lf",&energy);
      queue_fadd(qe,energy);
      }
    }
  /* convert queue to array - coordinates */
  if (qs->num) {
    if (g->geom)
      g->geom = gauss_geom_free(g->geom,g->n_geom,g->n_atoms);
    g->geom = gauss_geom_new(qs->num);
    if (!g->geom)
      msg_error("cannot allocate memory for array of structures",1);
    g->n_geom = 0;
    while (qs->num) {
      crd = queue_get(qs);
      g->geom[g->n_geom].coord = crd;
      /* convert to atomic units */
      for (i=0; i<g->n_atoms; i++)
        for (j=0; j<3; j++)
          g->geom[g->n_geom].coord[i][j] /= CONV_B_ANG;
      g->n_geom++;
      }
    }
  queue_free(qs);
  /* convert queue to array - energies */
  for (i=0; i<g->n_geom; i++) {
    if (qe->num)
      queue_fget(qe,&(g->geom[i].energy));
    else
      g->geom[i].energy = 0.0;
    }
  queue_free(qe);
  /* check if data was found */
  if (!g->n_geom) {
    if (c==GAUSS_CRD_INP)
      msg_error("cannot find input orientation coordinates",1);
    else if (c==GAUSS_CRD_STD)
      msg_error("cannot find standard orientation coordinates",1);
    else if (c==GAUSS_CRD_STD)
      msg_error("cannot find z-matrix orientation coordinates",1);
    else if (c==GAUSS_CRD_TRJ)
      msg_error("cannot find trajectory output coordinates",1);
    else if (c==GAUSS_CRD_ANY)
      msg_error("cannot find structure coordinates",1);
    else
      msg_error("invalid specification of coordinate orientation",1);
    }
  /* use last coordinate set */
  for (i=0; i<g->n_atoms; i++)
    for (j=0; j<3; j++)
      g->atom[i].coord[j] = g->geom[g->n_geom-1].coord[i][j];
  }

/* read coordinates from gaussian log file
 
   g - pointer to guassian data struct
   c - specify coordinate orientation
   f - name of the gaussian log file */
void gauss_log_read_coord(struct gauss_dat *g, short c, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_coord_f(g,c,file);
  file_close(file);
  }

/* read gaussian data from frequency analysis

   g - pointer to gaussian data struct
   f - pointer to open file stream */
void gauss_log_read_freq_f(struct gauss_dat *g, FILE *f) {
  char *line;
  unsigned i,j,k,l,n;
  short found = 0;
  rewind(f);
  /* transition dipoles */
  line = str_ffind_new(f,"VibFq2-Diag2:");
  if (line) {
    line = str_free(line);
    /* allocate memory */
    if (!g->freq)
      g->freq = gauss_freq_new(g->n_vib_modes);
    for (i=0; i<g->n_vib_modes; i++) {
      /* dipole derivative */
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of file while reading "
          "vibrational transition dipoles",1);
      if (!strstr(line,"Dipole derivative wrt mode"))
        msg_error("dipole derivative expected",1);
      g->freq[i].mu_d = vec_falloc(3);
      str_subst(line,'D','E');
      if (sscanf(line+32,"%lf%lf%lf",&(g->freq[i].mu_d[0]),
          &(g->freq[i].mu_d[1]),&(g->freq[i].mu_d[2]))!=3)
        msg_error("cannot read dipole derivative",1);
      line = str_free(line);
      /* contribution to vibrational polarizability */
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of file while reading "
          "contribution to vibrational polarizability",1);
      if (!strstr(line,"Vibrational polarizability contributions from mode"))
        msg_error("contribution to vibrational polarizability expected",1);
      g->freq[i].pol = vec_falloc(3);
      if (sscanf(line+55,"%lf%lf%lf",&(g->freq[i].pol[0]),
          &(g->freq[i].pol[1]),&(g->freq[i].pol[2]))!=3)
        msg_error("cannot read contribution to vibrational polarizability",1);
      line = str_free(line);
      /* electric dipole transition moment */
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of file while reading "
          "electric dipole transition moment",1);
      if (!strstr(line,"Electric dipole transition moment wrt mode"))
        msg_error("electric dipole transition moment expected",1);
      g->freq[i].mu_e = vec_falloc(3);
      str_subst(line,'D','E');
      if (sscanf(line+48,"%lf%lf%lf",&(g->freq[i].mu_e[0]),
          &(g->freq[i].mu_e[1]),&(g->freq[i].mu_e[2]))!=3)
        msg_error("cannot read electric dipole transition moment",1);
      line = str_free(line);
      /* magnetic dipole transition moment */
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of file while reading "
          "magnetic dipole transition moment",1);
      if (!strstr(line,"Magnetic dipole transition moment wrt mode"))
        msg_error("magnetic dipole transition moment expected",1);
      g->freq[i].mu_m = vec_falloc(3);
      str_subst(line,'D','E');
      if (sscanf(line+48,"%lf%lf%lf",&(g->freq[i].mu_m[0]),
          &(g->freq[i].mu_m[1]),&(g->freq[i].mu_m[2]))!=3)
        msg_error("cannot read magnetic dipole transition moment",1);
      line = str_free(line);
      }
    found = 1;
    }
  /* read frequencies and displacements */
  if (!found)
    rewind(f);
  line = str_ffind_new(f,"and normal coordinates:");
  if (line) {
    /* allocate memory */
    if (!g->freq)
      g->freq = gauss_freq_new(g->n_vib_modes);
    for (i=0; i<g->n_vib_modes; i++)
      if (!g->freq[i].displ)
        g->freq[i].displ = mat_falloc(g->n_atoms,3);
    /* data blocks */
    for (i=0; i<g->n_vib_modes/3; i++) {
      /* skip header */
      for (j=0; j<2; j++) {
        line = str_free(line);
        line = gauss_log_read_line("frequencies",f);
        }
      /* frequencies / intensities */
      n = ((i+1)<(g->n_vib_modes/3) ? 3 : 
           (g->n_vib_modes%3 ? g->n_vib_modes : 3));
      for (j=0, found=0; j<12; j++) {
        line = str_free(line);
        line = gauss_log_read_line("normal-mode data",f);
        /* Vibrational frequency */
        if (strstr(line," Frequencies --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].freq))!=1) 
              msg_error("cannot read frequency value",1);
          }
        /* Reduced mass */
        if (strstr(line," Red. masses --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].rmass))!=1) 
              msg_error("cannot read reduced-mass value f",1);
          }
        /* Force constant */
        if (strstr(line," Frc consts  --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].fconst))!=1) 
              msg_error("cannot read force-constant value f",1);
          }
        /* IR intensity */
        if (strstr(line," IR Inten    --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].ir))!=1) 
              msg_error("cannot read ir-intensity value f",1);
          }
        /* Raman activity */
        if (strstr(line," Raman Activ --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].raman))!=1) 
              msg_error("cannot read raman-activity value f",1);
          }
        /* Depolarization ratio (plane incident light) */
        if (strstr(line," Depolar (P) --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].depol_p))!=1) 
              msg_error("cannot read plane-depolarization-ratio value f",1);
          }
        /* Depolarization ratio (unpolarized incident light) */
        if (strstr(line," Depolar (U) --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].depol_u))!=1) 
              msg_error("cannot read unpolar-depolarization-ratio value f",1);
          }
        /* Dipole strength */
        if (strstr(line," Dip. str.   --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].strength_dip))!=1) 
              msg_error("cannot read dipole-strength value f",1);
          }
        /* Rotational strength */
        if (strstr(line," Rot. str.   --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].strength_rot))!=1) 
              msg_error("cannot read rotational-strength value f",1);
          }
        /* E-M angle */
        if (strstr(line," E-M angle   --")) {
          for (k=0; k<n; k++)
            if (sscanf(line+15+k*20,"%lf",&(g->freq[3*i+k].em_angle))!=1) 
              msg_error("cannot read e/m-angle value f",1);
          }
        /* Displacements */
        else if (strstr(line,"  Atom  AN")) {
          found = 1;
          break;
          }
        }
      if (!found)
        msg_error("cannot find vibrational displacements",1);
      /* displacements */
      for (j=0; j<g->n_atoms; j++) {
        line = str_free(line);
        line = gauss_log_read_line("vibrational displacements",f);
        for (k=0; k<n; k++)
          for (l=0; l<3; l++)
            if (sscanf(line+12+k*23+l*7,"%lf",&(g->freq[3*i+k].displ[j][l]))!=1)
              msg_error("cannot read vibrational displacements",1);
        }
      line = str_free(line);
      }
    }
  }

/* read gaussian data from frequency analysis

   g - pointer to gaussian data struct
   f - name of gaussian log file */
void gauss_log_read_freq(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_freq_f(g,file);
  file_close(file);
  }

/* convert list of centers to basis set data struct 
 
   g - pointer to gaussian data struct
   d - pointer to list with basis set data */
void gauss_log_read_basis_cb(struct gauss_dat *g, struct list *d) {
  unsigned i,id_c;
  struct basis *b = g->bs;
  struct basis_center *c;
  struct ldata *p = d->first;
  b->n_centers = d->num;
  b->center = basis_center_new(b->n_centers);
  id_c = 0;
  /* convert list */
  while (p) {
    c = (struct basis_center*)p->l_data;
    /* move data */
    b->center[id_c].type = c->type;
    b->center[id_c].n_shells = c->n_shells;
    for (i=0; i<3; i++)
      b->center[id_c].coord[i] = c->coord[i];
    b->center[id_c].shell = c->shell;
    p = p->l_next;
    id_c++;
    }
  /* clean memory */
  list_free(d,NULL);
  }

/* convert list of shells to list of centers
 
   g - pointer to gaussian data struct
   b - pointer to list with basis set data */
void gauss_log_read_basis_sc(struct gauss_dat *g, struct list *b) {
  unsigned i,id_c,id_s,n_shell;
  struct gauss_log_bs *d;
  struct basis_center *c;
  struct ldata *p1,*p2;
  struct list *store;
  n_shell = 0;
  p1 = p2 = b->first;
  store = list_alloc();
  id_c = ((struct gauss_log_bs*)p1->l_data)->at_id;
  /* create list of centers */
  while (p2) {
    d = (struct gauss_log_bs*)p2->l_data;
    n_shell++;
    if (!p2->l_next ||
       ((struct gauss_log_bs*)(p2->l_next->l_data))->at_id!=id_c) {
      /* collect shells from one center */
      c = basis_center_new(1);
      c->type = d->at_num;
      c->n_shells = n_shell;
      for (i=0; i<3; i++)
        c->coord[i] = d->coord[i];
      c->shell = basis_shell_new(c->n_shells);
      id_s = 0;
      /* move data from one list to another */
      while (p1!=p2->l_next) {
        d = (struct gauss_log_bs*)p1->l_data;
        c->shell[id_s].type = d->shell->type;
        c->shell[id_s].n_prim = d->shell->n_prim;
        c->shell[id_s].exp = d->shell->exp;
        c->shell[id_s].cf1 = d->shell->cf1;
        c->shell[id_s].cf2 = d->shell->cf2;
        p1 = p1->l_next;
        id_s++;
        }
      /* save to new list */
      list_add_end_p(store,c);
      id_c = d->at_id;
      n_shell = 0;
      }
    p2 = p2->l_next;
    }
  /* convert to basis set data struct */
  gauss_log_read_basis_cb(g,store);
  /* clean memory */
  list_free(b,NULL);
  }

/* read basis set from open gaussian log file
 
   g - pointer to gaussian data struct
   f - pointer to open file stream */
void gauss_log_read_basis_f(struct gauss_dat *g, FILE *f) {
  char *line,name[80],type[80],snum1[80],snum2[80];
  unsigned id,ip;
  struct gauss_log_bs *d;
  struct list *store;
  /* clean basis set struct */
  g->bs = basis_free(g->bs);
  g->bs = basis_new();
  /* find GFPrint format basis set */
  rewind(f);
  store = list_alloc();
  line = str_ffind_new(f,"AO basis set (Overlap normalization):");
  if (line) {
    line = str_free(line);
    line = gauss_log_read_line("basis set header",f);
    while (str_sub_bfind(line," Atom ")) {
      d = gauss_log_bs_new();
      if (sscanf(line,"%*s%s%*s%*u%s%u%*s%*s%*s%*s%lf%lf%lf",
        name,type,&(d->shell->n_prim),&(d->coord[0]),&(d->coord[1]),
        &(d->coord[2]))!=6)
        msg_error("invalid format of basis set, cannot read shell info",1);
      /* atom name and ID */
      id = 0;
      while (id<str_length(line) && isalpha(name[id]))
        id++;
      if (sscanf(name+id,"%u",&(d->at_id))!=1)
        msg_error("invalid format of atom label in basis set data",1);
      d->at_num = atom_num_pdb(name);
      /* primitive functions */
      d->shell->type = basis_shell_id(type,0);
      d->shell->exp = vec_falloc(d->shell->n_prim);
      d->shell->cf1 = vec_falloc(d->shell->n_prim);
      if (d->shell->type==BASIS_SHELL_SP)
        d->shell->cf2 = vec_falloc(d->shell->n_prim);
      for (ip=0; ip<d->shell->n_prim; ip++) {
        line = str_free(line);
        line = gauss_log_read_line("primitive basis set function",f);
        if (sscanf(line,"%s%s",snum1,snum2)!=2)
           msg_error("invalid format of primitive basis set function",1);
        d->shell->exp[ip] = str_fnum_d2e(snum1);
        d->shell->cf1[ip] = str_fnum_d2e(snum2);
        if (d->shell->type==BASIS_SHELL_SP) {
          if (sscanf(line+40,"%s",snum2)!=1)
            msg_error("invalid format of SP primitive basis set function",1);
          d->shell->cf2[ip] = str_fnum_d2e(snum2);
          }
        }
      list_add_end_p(store,d);
      line = str_free(line);
      line = gauss_log_read_line("basis set shell",f);
      }
    }
  /* convert data to basis set struct */
  if (store->num) {
    gauss_log_read_basis_sc(g,store);
    basis_def_parm(g->bs);
    }
  }

/* read basis set from gaussian log file
 
   g - pointer to gaussian data struct
   f - name of gaussian log file */
void gauss_log_read_basis(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_basis_f(g,file);
  file_close(file);
  }

/* read data from gaussian log file
 
   g - pointer to gaussian data struct
   f - name of the log file
   t - type of coordinates */
void gauss_log_read(struct gauss_dat *g, char *f, short t) {
  FILE *file;
  file = file_open(f,"r");
  gauss_log_read_spec_f(g,file);
  gauss_log_read_coord_f(g,t,file);
  gauss_log_read_basis_f(g,file);
  gauss_log_read_freq_f(g,file);
  gauss_log_read_states_f(g,file);
  gauss_log_read_nmr_f(g,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
