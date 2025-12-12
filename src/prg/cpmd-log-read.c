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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/units.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <mol/cell.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* read charge population from open cpmd log file
 
   d - pointer to cpmd data struct
   f - pointer to open cpmd file */
void cpmd_log_read_pop_f(struct cpmd_dat *d, FILE *f) {
  char *line;
  unsigned i;
  /* locate KS states record */
  rewind(f);
  line = str_ffind_new(f,"POPULATION ANALYSIS FROM PROJECTED");
  if (line) {
    /* read data */
    for (i=0; i<2; i++) {
      line = str_free(line);
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cpmd file while reading charges",1);
      }
    if (!strstr(line,"MULLIKEN") || !strstr(line,"LOWDIN"))
      msg_error("invalid format of population analysis",1);
    line = str_free(line);
    for (i=0; i<d->n_atoms; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cpmd file while reading charges",1);
      if (sscanf(line,"%*d%*s%lf%lf%lf",&(d->atom[i].charge[CPMD_CHARGE_MUL]),
        &(d->atom[i].charge[CPMD_CHARGE_LOW]),&(d->atom[i].valence))!=3)
        msg_error("invalid format of population analysis",1);
      line = str_free(line);
      }
    }
  line = str_free(line);
  }

/* read charge population from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_pop(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_pop_f(d,file);
  file_close(file);
  }

/* read atom specification from cpmd output file
 
   d - pointer to cpmd data struct
   f - pointer to open cpmd file */
void cpmd_log_read_atoms(struct cpmd_dat *d, FILE *f) {
  char *line,sym[80];
  unsigned i,id = 0,at[120];
  struct cpmd_atom *a;
  struct queue *q;
  /* initalization */
  vec_uset(at,0,120);
  /* skip header */
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of cpmd file while reading atoms",1);
  line = str_free(line);
  /* read atom types and coordinates */
  q = queue_alloc();
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (line[1]=='*')
      break;
    a = cpmd_atom_new(1);
    if (sscanf(line,"%*d%s%lf%lf%lf",sym,&(a->coord[0]),
      &(a->coord[1]),&(a->coord[2]))!=4) {
      cpmd_atom_free(a,1);
      msg_error("invalid format of atom specification in cpmd file",1);
      }
    a->t_id = atom_num(sym);
    vec_fset(a->charge,0.0,3);
    a->valence = 0.0;
    at[a->t_id] = 1;
    queue_add(q,a);
    }
  /* save data to cpmd struct */
  d->n_atoms = q->num;
  d->atom = cpmd_atom_new(d->n_atoms);
  while (q->num) {
    a = queue_get(q);
    for (i=0; i<3; i++)
      d->atom[id].coord[i] = a->coord[i]*CONV_B_ANG;
    d->atom[id].t_id = a->t_id;
    cpmd_atom_free(a,1);
    id++;
    }
  /* set number of types */
  cpmd_dat_set_types(d,at);
  /* clean memory */
  queue_free(q);
  }

/* read kohn-sham eigenvalues from cpmd output file - one state
 
   v - array for energy state data
   m - number of states
   s - spin polarized calculation 
   f - pointer to open file */
void cpmd_log_read_states_one(struct cpmd_state *v, unsigned m,
  short s, FILE *f) {
  char *line,num1[256],num2[256];
  unsigned i,n = m/2;
  if (s) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cpmd output file while reading"
        " LSD states header ",1);
    line = str_free(line);
    }
  for (i=0; i<n; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cpmd output file while reading states",1);
    if (sscanf(line,"%*d%lf%lf%*d%lf%lf",&(v[2*i].energy),&(v[2*i].occup),
      &(v[2*i+1].energy),&(v[2*i+1].occup))!=4) {
      if (strstr(line,"NC")) {
        if (sscanf(line,"%*d%s%lf%*d%s%lf",num1,&(v[2*i].occup), 
          num2,&(v[2*i+1].occup))!=4)
          msg_error("problem with reading non-converged states",1);
        if (sscanf(num1,"%lfNC",&(v[2*i].energy))!=1)
          msg_error("problem with reading non-converged states - column 1",1);
        if (sscanf(num2,"%lfNC",&(v[2*i+1].energy))!=1)
          msg_error("problem with reading non-converged states - column 2",1);
        }
      else
        msg_error("invalid format of states in cpmd output file",1);
      }
    line = str_free(line);
    }
  if (m%2) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cpmd file while reading last state",1);
    if (sscanf(line,"%*d%lf%lf",&(v[m-1].energy),&(v[m-1].occup))!=2) {
      if (strstr(line,"NC")) {
        if (sscanf(line,"%*d%s%lf",num1,&(v[m-1].occup))!=2)
          msg_error("problem with reading non-converged states",1);
        if (sscanf(num1,"%lfNC",&(v[m-1].energy))!=1)
          msg_error("problem with reading non-converged states",1);
        }
      else
        msg_error("invalid format of last state in cpmd output file",1);
      }
    line = str_free(line);
    }
  for (i=0; i<m; i++) {
    v[i].compln = 0.0;
    v[i].ao_pj = NULL;
    }
  }

/* read state data from open cpmd log file
 
   d - pointer to cpmd data struct
   f - pointer to open file */
void cpmd_log_read_states_f(struct cpmd_dat *d, FILE *f) {
  char *line;
  short kpoints = 0;
  unsigned i;
  fpos_t fpos;
  /* locate KS states record */
  rewind(f);
  line = str_ffind_new(f,"EIGENVALUES(EV) AND OCCUPATION:");
  if (line) {
    line = str_free(line);
    /* read data */
    fgetpos(f,&fpos);
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cpmd file while reading states",1);
    if (strstr(line,"K POINT:")) {
      kpoints = 1;
      for (line=str_free(line),line=str_read_line_new(f); line;
           line=str_free(line),line=str_read_line_new(f)) {
        if (strstr(line,"CHEMICAL POTENTIAL =")) {
          line = str_free(line);
          break;
          }
        if (strstr(line,"K POINT:"))
          kpoints++;
        }
      }
    fsetpos(f,&fpos);
    /* K points */
    if (kpoints) {
      d->n_kpoints = kpoints;
      if (d->lsd_calc) 
        msg_error("LSD calculation of K-points is not supported yet",1);
      else {
        d->state_a = cpmd_state_new(d->n_states*d->n_kpoints);
        for (i=0; i<d->n_kpoints; i++) {
          line = str_read_line_new(f);
          if (!line)
            msg_error("unexpected end of file while reading"
              " k point energies",1);
          line = str_free(line);
          cpmd_log_read_states_one(d->state_a+d->n_states*i,d->n_states,0,f);
          }
        }
      }
    /* single KS */
    else {
      if (d->lsd_calc) {
        d->state_a = cpmd_state_new(d->n_alpha_states);
        cpmd_log_read_states_one(d->state_a,d->n_alpha_states,1,f);
        d->state_b = cpmd_state_new(d->n_beta_states);
        cpmd_log_read_states_one(d->state_b,d->n_beta_states,1,f);
        }
      else {
        d->state_a = cpmd_state_new(d->n_states);
        cpmd_log_read_states_one(d->state_a,d->n_states,0,f);
        }
      }
    }
  }

/* read state data from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_states(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_states_f(d,file);
  file_close(file);
  }

/* read projection of one set of states to atomic basis set
 
   v       - array for projection data
   n_state - number of states
   n_orb   - number of orbitals
   s       - spin polarized calculation 
   f       - pointer to open cpmd file */
void cpmd_log_read_proj_one(struct cpmd_state *v, unsigned n_state,
  unsigned n_orb, short s, FILE *f) {
  char msg[128] = "unexpected end of cpmd file while reading projection";
  unsigned i,ia,ib,is,nb,ns,off;
  char *line;
  if (s) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cpmd output file while reading"
        " LSD projection of states header ",1);
    line = str_free(line);
    }
  /* read projection coefficients */
  line = str_read_line_new(f);
  if (!line)
    msg_error(msg,1);
  line = str_free(line);
  nb = (n_state%8 ? n_state/8+1 : n_state/8);
  for (ib=0; ib<nb; ib++) {
    /* skip header */
    for (i=0; i<2; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error(msg,1);
      line = str_free(line);
      }
    /* projection coefficients */
    ns = (ib==nb-1 ? (n_state%8 ? n_state%8 : 8) : 8);
    ia = 0;
    while (ia<n_orb) {
      line = str_read_line_new(f);
      if (!line)
        msg_error(msg,1);
      if (str_length(line)>10) {
        off = 14;
        if (strstr(line,"COMPLETNESS")) {
          for (is=0; is<ns; is++)
            if (sscanf(line+off+is*8,"%lf",&(v[ib*8+is].compln))!=1)
              msg_error("invalid format of projection data - completness",1);
          }
        else if (strstr(line,"OCCUPATION")) {
          for (is=0; is<ns; is++)
            if (sscanf(line+off+is*8,"%lf",&(v[ib*8+is].occup))!=1)
              msg_error("invalid format of projection data - occupation",1);
          }
        else {
          for (is=0; is<ns; is++)
            if (sscanf(line+off+is*8,"%lf",&(v[ib*8+is].ao_pj[ia]))!=1)
              msg_error("invalid format of projection data - coefficients",1);
          ia++;
          }
        }
      line = str_free(line);
      }
    }
  }

/* read projection of states from open cpmd log file
 
   d - pointer to cpmd data struct
   f - pointer to open cpmd file */
void cpmd_log_read_proj_f(struct cpmd_dat *d, FILE *f) {
  char *line;
  unsigned i;
  /* locate KS states record */
  rewind(f);
  line = str_ffind_new(f,"WAVEFUNCTIONS IN ATOMIC ORBITAL BASIS");
  if (line) {
    line = str_free(line);
    /* spin polarized system */
    if (d->lsd_calc) {
      if (!d->state_a) {
        d->state_a = cpmd_state_new(d->n_alpha_states);
        for (i=0; i<d->n_alpha_states; i++) {
          d->state_a[i].energy = 0.0;
          d->state_a[i].compln = 0.0;
          d->state_a[i].occup = 0.0;
          }
        }
      for (i=0; i<d->n_alpha_states; i++)
        d->state_a[i].ao_pj = vec_falloc(d->n_orbitals);
      cpmd_log_read_proj_one(d->state_a,d->n_alpha_states,
        d->n_orbitals,1,f);
      if (!d->state_b) {
        d->state_b = cpmd_state_new(d->n_beta_states);
        for (i=0; i<d->n_beta_states; i++) {
          d->state_b[i].energy = 0.0;
          d->state_b[i].compln = 0.0;
          d->state_b[i].occup = 0.0;
          }
        }
      for (i=0; i<d->n_beta_states; i++)
        d->state_b[i].ao_pj = vec_falloc(d->n_orbitals);
      cpmd_log_read_proj_one(d->state_b,d->n_beta_states,
        d->n_orbitals,1,f);
      }
    /* closed shell system */
    else {
      if (!d->state_a) {
        d->state_a = cpmd_state_new(d->n_states);
        for (i=0; i<d->n_states; i++) {
          d->state_a[i].energy = 0.0;
          d->state_a[i].compln = 0.0;
          d->state_a[i].occup = 0.0;
          }
        }
      for (i=0; i<d->n_states; i++)
        d->state_a[i].ao_pj = vec_falloc(d->n_orbitals);
      cpmd_log_read_proj_one(d->state_a,d->n_states,
        d->n_orbitals,0,f);
      }
    }
  }

/* read projection of states from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_proj(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_proj_f(d,file);
  file_close(file);
  }

/* convert spin multiplicity keyword to numerical value
 
   s - the keyword */
short cpmd_log_read_mult(char *s) {
  if (str_compare(s,"SINGLET")) return(1);
  if (str_compare(s,"DOUBLET")) return(2);
  if (str_compare(s,"TRIPLET")) return(3);
  msg_error_f("unknown spin multiplicity \"%s\"",1,s);
  return(0);
  }

/* read molecule specification from open cpmd log file
 
   d - pointer ot cpmd data struct
   f - pointer ot open cpmd log file */
void cpmd_log_read_spec_f(struct cpmd_dat *d, FILE *f) {
  char *line,str[1024];
  double lv[3][3];
  unsigned i,j;
  fpos_t fpos;
  /* atomic coordinates */
  line = str_ffind_new(f,"**** ATOMS ****");
  if (!line)
    msg_error("cannot found atom specification in cpmd output",1);
  line = str_free(line);
  cpmd_log_read_atoms(d,f);
  /* k points */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (strstr(line,"SPECIAL K-POINTS GENERATION")) {
      for (i=0; i<3; i++) {
        line = str_free(line);
        line = str_read_line_new(f);
        }
      if (sscanf(line+40,"%u%u%u",&(d->kp_mesh[0]),
           &(d->kp_mesh[1]),&(d->kp_mesh[2]))!=3)
        msg_error("cannot read k-point mesh dimensions",1);
      for (i=0; i<5; i++) {
        line = str_free(line);
        line = str_read_line_new(f);
        }
      if (sscanf(line+56,"%u",&(d->n_kpoints))!=1)
        msg_error("cannot read number of special k-points",1);
      line = str_free(line);
      line = str_read_line_new(f);
      d->kpoint = cpmd_kpoint_new(d->n_kpoints);
      for (i=0; i<d->n_kpoints; i++) {
        line = str_free(line);
        line = str_read_line_new(f);
        if (sscanf(line,"%d%lf%lf%lf%lf",&j,
            &(d->kpoint[i].coord[0]),&(d->kpoint[i].coord[1]),
            &(d->kpoint[i].coord[2]),&(d->kpoint[i].weight))!=5)
          msg_error_f("cannot read k-point #%d specification",1,i+1);
        }
      }
    if (str_sub_bfind(line," NUMBER OF STATES:"))
      break;
    }
  /* number of states */
  if (!line)
    msg_error("cannot found number of states in cpmd output",1);
  if (sscanf(line+18,"%u",&(d->n_states))!=1)
    msg_error("cannot read number of states in cpmd output",1);
  line = str_free(line);
  /* number of electrons */
  line = str_ffind_new(f,"NUMBER OF ELECTRONS:");
  if (!line)
    msg_error("cannot found number of electrons in cpmd output",1);
  if (sscanf(line+21,"%u",&(d->n_electrons))!=1)
    msg_error("cannot read number of electrons in cpmd output",1);
  line = str_free(line);
  /* system charge */
  line = str_ffind_new(f,"CHARGE:");
  if (!line)
    msg_error("cannot found charge in cpmd output",1);
  if (sscanf(line+8,"%lf",&(d->charge))!=1)
    msg_error("cannot read charge in cpmd output",1);
  line = str_free(line);
  /* open vs. closed shell */
  fgetpos(f,&fpos);
  line = str_ffind_new(f,"SPIN MULTIPLICITY");
  if (line) {
    line = str_free(line);
    d->lsd_calc = 1;
    if (sscanf(line+20,"%s",str)!=1)
      msg_error("cannot read spin multiplicity in cpmd output",1);
    d->multiplicity = cpmd_log_read_mult(str);
    /* number of alpha states */
    line = str_ffind_new(f,"NUMBER OF ALPHA STATES");
    if (!line)
      msg_error("cannot found number of alpha states in cpmd output",1);
    if (sscanf(line+25,"%u",&(d->n_alpha_states))!=1)
      msg_error("cannot read number of alpha states in cpmd output",1);
    line = str_free(line);
    /* number of beta states */
    line = str_ffind_new(f,"NUMBER OF BETA STATES");
    if (!line)
      msg_error("cannot found number of beta states in cpmd output",1);
    if (sscanf(line+25,"%u",&(d->n_beta_states))!=1)
      msg_error("cannot read number of beta states in cpmd output",1);
    line = str_free(line);
    }
  else {
    fsetpos(f,&fpos);
    d->lsd_calc = 0;
    }
  /* lattice vectors */
  line = str_ffind_new(f,"LATTICE VECTOR A1");
  if (!line)
    msg_error("cannot found lattice vector a1 in cpmd output",1);
  if (sscanf(line+25,"%lf%lf%lf",&(lv[0][0]),&(lv[0][1]),&(lv[0][2]))!=3)
    msg_error("invalid format of lattice vector a1 in cpmd output",1);
  line = str_free(line);
  line = str_ffind_new(f,"LATTICE VECTOR A2");
  if (!line)
    msg_error("cannot found lattice vector a2 in cpmd output",1);
  if (sscanf(line+25,"%lf%lf%lf",&(lv[1][0]),&(lv[1][1]),&(lv[1][2]))!=3)
    msg_error("invalid format of lattice vector a2 in cpmd output",1);
  line = str_free(line);
  line = str_ffind_new(f,"LATTICE VECTOR A3");
  if (!line)
    msg_error("cannot found lattice vector a3 in cpmd output",1);
  if (sscanf(line+25,"%lf%lf%lf",&(lv[2][0]),&(lv[2][1]),&(lv[2][2]))!=3)
    msg_error("invalid format of lattice vector a3 in cpmd output",1);
  line = str_free(line);
  /* lattice parameters */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      lv[i][j]*=CONV_B_ANG;
  cell_set_vec(d->cell,lv[0],lv[1],lv[2]);
  }

/* read molecule specification from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_spec(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_spec_f(d,file);
  file_close(file);
  }

/* read atomic basis for from cpmd log file
 
   d - pointer to cpmd data struct
   f - pointer to open file */
void cpmd_log_read_basis_f(struct cpmd_dat *d, FILE *f) {
  char *line,sym[10];
  unsigned i,n,type_id = 0;
  struct queue *q;
  /* locate basis set record */
  rewind(f);
  line = str_ffind_new(f,"GENERATE ATOMIC BASIS SET");
  if (line) {
    line = str_free(line);
    q = queue_alloc();
    /* read data to queue */
    for (line=str_read_line_new(f); line;
         line=str_free(line),line=str_read_line_new(f)) {
      if (str_length(line)<10)
        break;
      if (strstr(line,"OCCUPATION")) {
        if (strstr(line,"ALPHA")) {
          if (sscanf(line,"%s",sym)!=1)
            msg_error("cannot read orbital type in basis set specification",1);
          }
        else if (strstr(line,"L VALUE")) {
          if (sscanf(line+16,"%s",sym)!=1)
            msg_error("cannot read orbital type in basis set specification",1);
          }
        else
          msg_error("unknown basis set specification",1);
        queue_siadd(q,cpmd_bs_id(sym));
        }
      else if (strstr(line," ORBITALS")) {
        if (sscanf(line,"%s",sym)!=1)
          msg_error("cannot read atom type in cpmd basis set specification",1);
        cpmd_bs_set(d,q,type_id);
        n = atom_num(sym);
        for (i=0; i<d->n_types; i++)
          if (d->type[i].num==n)
            break;
        type_id = i;
        }
      else
        msg_error("format of basis set in cpmd file not recognized",1);
      }
    /* convert queue to cpmd basis set struct array */
    cpmd_bs_set(d,q,type_id);
    queue_free(q);
    /* set total number of orbitals */
    for (i=0; i<d->n_atoms; i++)
      d->n_orbitals += (d->type[d->atom[i].t_id].n_orbs);
    }
  }

/* read atomic basis set from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_basis(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_basis_f(d,file);
  file_close(file);
  }

/* read system structure from open cpmd log file
 
   d - pointer to cpmd data struct
   f - pointer to open file */
void cpmd_log_read_geom_f(struct cpmd_dat *d, FILE *f) {
  char *line,word[80],sym[80];
  unsigned i,j,t;
  double **crd;
  long id;
  struct queue *q;
  /* locate geometry optimization */
  rewind(f);
  line = str_ffind_new(f,"      GEOMETRY OPTIMIZATION      ");
  if (line) {
    line = str_free(line);
    /* read structures to queue */
    q = queue_alloc();
    for (line=str_ffind_new(f,"COORDINATES            GRADIENTS (-FORCES)");
         line; line=str_free(line),
         line=str_ffind_new(f,"COORDINATES            GRADIENTS (-FORCES)")) {
      crd = mat_falloc(d->n_atoms+1,3);
      /* coordinates */
      for (i=0; i<d->n_atoms; i++) {
        line = str_read_line_new(f);
        if (!line)
          break;
        str_sub_copy(line,word,7,15);
        if (sscanf(word,"%lf",&(crd[i][0]))!=1)
          break;
        str_sub_copy(line,word,15,23);
        if (sscanf(word,"%lf",&(crd[i][1]))!=1)
          break;
        str_sub_copy(line,word,23,31);
        if (sscanf(word,"%lf",&(crd[i][2]))!=1)
          break;
        line = str_free(line);
        }
      /* energy */
      for (j=0; j<3; j++) {
        line = str_read_line_new(f);
        if (!line)
          break;
        line = str_free(line);
        }
      if (!line)
        break;
      if (sscanf(line+46,"%lf",&(crd[d->n_atoms][0]))!=1)
        break;
      queue_add(q,crd);
      }
    /* convert queue to structure array */
    d->n_geoms = q->num;
    d->g_opt = cpmd_geom_new(d->n_geoms,d->n_atoms);
    id = 0;
    while (q->num) {
      crd = queue_get(q);
      for (i=0; i<d->n_atoms; i++)
        for (j=0; j<3; j++)
          d->g_opt[id].coord[i][j] = (crd[i][j]*CONV_B_ANG);
      d->g_opt[id].energy = crd[d->n_atoms][0];
      mat_ffree(crd,d->n_atoms+1);
      id++;
      }
    queue_free(q);
    }
  /* locate final structure */
  rewind(f);
  line = str_ffind_new(f,"FINAL RESULTS");
  if (line) {
    /* header */
    for (i=0; i<4; i++) {
      line = str_free(line);
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of file while reading final geometry",1);
      }
    /* geometry optimization format */
    if (strstr(line,"COORDINATES            GRADIENTS (-FORCES)")) {
      /* allocate memory for coordinates */
      d->g_final = cpmd_geom_new(1,d->n_atoms);
      /* read data to queue */
      for (i=0; i<d->n_atoms; i++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of file while reading final geometry",1);
        str_sub_copy(line,word,7,15);
        if (sscanf(word,"%lf",&(d->g_final->coord[i][0]))!=1)
          msg_error("invalid format of final geometry coordinates",1);
        str_sub_copy(line,word,15,23);
        if (sscanf(word,"%lf",&(d->g_final->coord[i][1]))!=1)
          msg_error("invalid format of final geometry coordinates",1);
        str_sub_copy(line,word,23,31);
        if (sscanf(word,"%lf",&(d->g_final->coord[i][2]))!=1)
          msg_error("invalid format of final geometry coordinates",1);
        for (j=0; j<3; j++)
          d->g_final->coord[i][j]*=CONV_B_ANG;
        line = str_free(line);
        }
      }
    /* single point format */
    else {
      line = str_ffind_new(f,"ATOMIC COORDINATES");
      if (line) {
        line = str_free(line);
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of file while reading final geometry",1);
        line = str_free(line);
        /* allocate memory for coordinates */
        d->g_final = cpmd_geom_new(1,d->n_atoms);
        /* read data to queue */
        for (i=0; i<d->n_atoms; i++) {
          line = str_read_line_new(f);
          if (!line)
            msg_error("unexpected end of file while reading final geometry",1);
          if (sscanf(line+10,"%s%lf%lf%lf",sym,
            &(d->g_final->coord[i][0]),
            &(d->g_final->coord[i][1]),
            &(d->g_final->coord[i][2]))!=4)
            msg_error("invalid format of final structure coordinates",1);
          t = atom_num(sym);
          if (d->type[d->atom[i].t_id].num!=t)
            msg_error("incompatible atom types in the final structure",1);
          for (j=0; j<3; j++)
            d->g_final->coord[i][j]*=CONV_B_ANG;
          line = str_free(line);
          }
        }
      /* not found */
      else
        msg_error("cannot locate final structure coordinates",1);
      }
    }
  }

/* read system structure from cpmd log file
 
   d - pointer ot cpmd data struct
   f - name of the log file */
void cpmd_log_read_geom(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_geom_f(d,file);
  file_close(file);
  }

/* read data from cpmd output file
 
   d - pointer to cpmd data struct
   f - name of the file */
void cpmd_log_read(struct cpmd_dat *d, char *f) {
  FILE *file;
  file = file_open(f,"r");
  cpmd_log_read_spec_f(d,file);
  cpmd_log_read_basis_f(d,file);
  cpmd_log_read_states_f(d,file);
  cpmd_log_read_proj_f(d,file);
  cpmd_log_read_pop_f(d,file);
  cpmd_log_read_geom_f(d,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
