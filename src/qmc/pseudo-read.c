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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <cmn/file.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "qmc/pseudo.h"

/* -------------------------------------------------------------------------- */

/* assign value from file to internal variable
 
   s - line from from pseudopotential file
   d - pointer to internal variable
   t - data type */
void pseudo_read_value(char *s, void *d, short t) {
  unsigned m = 0;
  switch (t) {
    case TYPE_UINT:   m = sscanf(s,"%u", (unsigned*)d); break;
    case TYPE_DOUBLE: m = sscanf(s,"%lf",(double*)d);   break;
    }
  if (m!=1)
    msg_error_f("unreadable value in pseudopotential file: \"%s\"",1,s);
  }

/* read pseudopotential data from file - atom section

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_atom(struct pseudo *p, FILE *f) {
  char *line;
  line = str_ffind_new(f,"&ATOM");
  if (!line)
    msg_error("cannot found atom section in pseudopotential file",1); 
  line = str_free(line);
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (str_sub_bfind(line,"&END"))
      break;
    else if (strstr(line,"Z  ="))
      pseudo_read_value(line+5,&(p->atom_num),TYPE_UINT);
    else if (strstr(line,"ZV ="))
      pseudo_read_value(line+5,&(p->val_num),TYPE_UINT);
    else if (strstr(line,"XC =")) {
      pseudo_read_value(line+5,&(p->xc_num),TYPE_UINT);
      pseudo_read_value(line+12,&(p->xc_slater),TYPE_DOUBLE);
      }
    else if (strstr(line,"TYPE =")) {
      if (strstr(line+7,"NORMCONSERVING NUMERIC"))
        p->pp_type = PSEUDO_NC_NUM;
      else
        msg_error("unknown pseudopotential type",1);
      }
    }
  }

/* read pseudopotential data from file - electron configuration

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_econf(struct pseudo *p, FILE *f) {
  char *line,sym[10];
  unsigned i;
  for (i=0; i<(p->state_core+p->state_val); i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of pseudopotential file",1);
    if (sscanf(line+30,"%d%s%lf",&(p->el_conf_n[i]),sym,
      &(p->el_conf_occ[i]))!=3)
      msg_error("invalid format of electron configuration",1);
    p->el_conf_l[i] = pseudo_shell_id(sym);
    }
  }

/* read pseudopotential data from file - exchange-correlation functional

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_func(struct pseudo *p, FILE *f) {
  char *line;
  double num;
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading slater exchange",1);
  pseudo_read_value(line+27,&(num),TYPE_DOUBLE);
  if (fabs(p->xc_slater-num)>0.001)
    msg_error("inconsistent slater exchange in pseudopotential file",1);
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading lda correlation",1);
  p->xc_lda_c = pseudo_func_id(line+27);
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading gc exchange",1);
  p->xc_gc_e = pseudo_func_id(line+27);
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading gc correlation",1);
  p->xc_gc_c = pseudo_func_id(line+27);
  }

/* read trouiller-martins pseudopotential data 

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_mt_pp(struct pseudo *p, FILE *f) {
  char *line,sym[10];
  double pp_r[100],pp_e[100];
  unsigned n=0,pp_n[100];
  short pp_l[100];
  fpos_t fpos;
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading mt pp data",1);
  line = str_free(line);
  fgetpos(f,&fpos);
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (sscanf(line+5,"%d%s%lf%lf",&(pp_n[n]),sym,&(pp_r[n]),&(pp_e[n]))!=4)
      break;
    pp_l[n] = pseudo_shell_id(sym);
    n++;
    }
  fsetpos(f,&fpos);
  p->mt_pp_n = vec_ucopy_new(pp_n,n);
  p->mt_pp_l = vec_sicopy_new(pp_l,n);
  p->mt_pp_r = vec_fcopy_new(pp_r,n);
  p->mt_pp_e = vec_fcopy_new(pp_e,n);
  p->mt_num = n;
  }

/* read pseudopotential data from file - info section

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_info(struct pseudo *p, FILE *f) {
  char *line,text[1024];
  unsigned num;
  line = str_ffind_new(f,"&INFO");
  if (!line)
    msg_error("cannot found info section in pseudopotential file",1); 
  line = str_free(line);
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (str_sub_bfind(line,"&END"))
      break;
    else if (strstr(line,"Pseudopotential Report")) {
      str_sub_copy(line,text,32,63);
      str_trim(text);
      p->time_stamp = str_copy_new(text);
      }
    else if (strstr(line,"Atomic Symbol")) {
      str_sub_copy(line,text,41,63);
      str_trim(text);
      if (p->atom_num!=atom_num(text))
        msg_error("inconsistent atomic symbol in info section",1);
      }
    else if (strstr(line,"Atomic Number")) {
      pseudo_read_value(line+40,&num,TYPE_UINT);
      if (p->atom_num!=num)
        msg_error("inconsistent atomic number in info section",1);
      }
    else if (strstr(line,"Number of core states"))
      pseudo_read_value(line+40,&(p->state_core),TYPE_UINT);
    else if (strstr(line,"Number of valence states")) {
      pseudo_read_value(line+40,&(p->state_val),TYPE_UINT);
      p->el_conf_n = vec_ualloc(p->state_core+p->state_val);
      p->el_conf_l = vec_ualloc(p->state_core+p->state_val);
      p->el_conf_occ = vec_falloc(p->state_core+p->state_val);
      }
    else if (strstr(line,"Exchange-Correlation Functional"))
      pseudo_read_func(p,f);
    else if (strstr(line,"Electron Configuration"))
      pseudo_read_econf(p,f);
    else if (strstr(line,"Full Potential Total Energy"))
      pseudo_read_value(line+36,&(p->en_pot),TYPE_DOUBLE);
    else if (strstr(line,"Trouiller-Martins normconserving PP"))
      pseudo_read_mt_pp(p,f);
    else if (strstr(line,"Number of Mesh Points"))
      pseudo_read_value(line+30,&(p->mesh_num),TYPE_UINT);
    else if (strstr(line,"Pseudoatom Total Energy"))
      pseudo_read_value(line+32,&(p->en_tot),TYPE_DOUBLE);
    }
  }

/* read pseudopotential data from file - potential section

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_potential(struct pseudo *p, FILE *f) {
  char *line,**sw;
  unsigned i,j,num = 0,sn;
  line = str_ffind_new(f,"&POTENTIAL");
  if (!line)
    msg_error("cannot found potential section in pseudopotential file",1); 
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading potential mesh",1);
  if (sscanf(line,"%d",&num)!=1)
    msg_error("invalid format of potential mesh data",1);
  if (p->mesh_num!=num)
    msg_error("inconsistent mesh number for potential",1);
  p->dat_pot = mat_falloc(p->mesh_num,p->mt_num+1);
  line = str_free(line);
  for (i=0; i<p->mesh_num; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of file while reading potential data",1);
    sw = str_split(line,' ',&sn);
    if (sn<(p->mt_num+1))
      msg_error("invalid format of potential data",1);
    if (sscanf(sw[0],"%lf",&(p->dat_pot[i][0]))!=1)
      msg_error_f("unrecognizable distance value \"%s\"",1,sw[0]);
    for (j=1; j<=p->mt_num; j++)
      if (sscanf(sw[j],"%lf",&(p->dat_pot[i][j]))!=1)
        msg_error_f("unrecognizable potential value \"%s\"",1,sw[j]);
    vec_sfree(sw,sn);
    str_free(line);
    }
  }

/* read pseudopotential data from file - wavefunction section

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_wavefunction(struct pseudo *p, FILE *f) {
  char *line,**sw;
  unsigned i,j,num = 0,sn;
  line = str_ffind_new(f,"&WAVEFUNCTION");
  if (!line)
    msg_error("cannot found wavefunction section in pseudopotential file",1); 
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading wavefunction mesh",1);
  if (sscanf(line,"%d",&num)!=1)
    msg_error("invalid format of wavefunction mesh data",1);
  if (p->mesh_num!=num)
    msg_error("inconsistent mesh number for wavefunction",1);
  p->dat_wfce = mat_falloc(p->mesh_num,p->mt_num+1);
  line = str_free(line);
 for (i=0; i<p->mesh_num; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of file while reading wavefunction data",1);
    sw = str_split(line,' ',&sn);
    if (sn<(p->mt_num+1))
      msg_error("invalid format of wavefunction data",1);
    if (sscanf(sw[0],"%lf",&(p->dat_wfce[i][0]))!=1)
      msg_error_f("unrecognizable distance value \"%s\"",1,sw[0]);
    for (j=1; j<=p->mt_num; j++)
      if (sscanf(sw[j],"%lf",&(p->dat_wfce[i][j]))!=1)
        msg_error_f("unrecognizable wavefunction value \"%s\"",1,sw[j]);
    vec_sfree(sw,sn);
    str_free(line);
    }
  }

/* read pseudopotential data from file - atom density section

   p - pointer to pseudopotential data struct
   f - pointer to open input file */
void pseudo_read_density(struct pseudo *p, FILE *f) {
  char *line,**sw;
  unsigned i,num = 0,sn;
  line = str_ffind_new(f,"&ATDENS");
  if (!line)
    return;
  line = str_free(line);
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of file while reading density mesh",1);
  if (sscanf(line,"%d",&num)!=1)
    msg_error("invalid format of density mesh data",1);
  if (p->mesh_num!=num)
    msg_error("inconsistent mesh number for density",1);
  p->dat_den = mat_falloc(p->mesh_num,2);
  line = str_free(line);
  for (i=0; i<p->mesh_num; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of file while reading density data",1);
    sw = str_split(line,' ',&sn);
    if (sn!=2)
      msg_error("invalid format of density data",1);
    if (sscanf(sw[0],"%lf",&(p->dat_den[i][0]))!=1)
      msg_error_f("unrecognizable distance value \"%s\"",1,sw[0]);
    if (sscanf(sw[1],"%lf",&(p->dat_wfce[i][1]))!=1)
      msg_error_f("unrecognizable density value \"%s\"",1,sw[1]);
    vec_sfree(sw,sn);
    str_free(line);
    }
  }

/* read pseudopotential data from file

   p - pointer to pseudopotential data struct
   f - name of the input file */
void pseudo_read(struct pseudo *p, char *f) {
  FILE *file;
  file = file_open(f,"r");
  pseudo_read_atom(p,file);
  pseudo_read_info(p,file);
  pseudo_read_potential(p,file);
  pseudo_read_wavefunction(p,file);
  pseudo_read_density(p,file);
  file_close(file); 
  }

/* -------------------------------------------------------------------------- */
