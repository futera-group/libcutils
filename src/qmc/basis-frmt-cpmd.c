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
#include <stdlib.h>
#include <cmn/file.h>
#include <cmn/list.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* read one non empty line from basis set file in CPMD format
 
   line - pointer to previously allocated line
   desc - description of read data
   file - pointer to open file stream */
char* basis_cpmd_read_line(char *line, char *desc, FILE *file) {
  line = str_free(line);
  line = str_read_line_new(file);
  if (!line && desc)
    msg_error_f("unexpected end of file while reading %s",1,desc);
  return(line);
  }

/* read basis set in CPMD format from file
 
   b - pointer to basis set data struct
   f - name of the file */
void basis_cpmd_read(struct basis *b, char *f) {
  char *line,name[256],type[256];
  unsigned i,j,k,id_c,id_s,id_t,id_p0,id_p1,n,n_types,n_shells,n_prims;
  double *de,*dc;
  FILE *file = NULL;
  struct list *st_c,*st_s;
  struct basis_center *c;
  struct basis_shell *s;
  struct ldata *p;
  if (f && f[0])
    file = file_open(f,"r");
  /* read data from file */
  line = NULL;
  st_c = list_alloc();
  for (line=basis_cpmd_read_line(line,NULL,file); line;
       line=basis_cpmd_read_line(line,NULL,file)) {
    /* new center */
    c = basis_center_new(1);
    if (sscanf(line,"%s",name)!=1)
      msg_error("invalid format of basis set header in CPMD basis set file",1);
    c->type = atom_num(name);
    /* shell types */
    line = basis_cpmd_read_line(line,"number of shell types",file);
    if (sscanf(line,"%u",&(n_types))!=1)
      msg_error("cannot read number of shell types in CPMD basis set file",1);
    /* read shells */
    st_s = list_alloc();
    for (id_t=0; id_t<n_types; id_t++) {
      /* shell type */
      line = basis_cpmd_read_line(line,"shell type",file);
      if (sscanf(line,"%s",type)!=1)
        msg_error("invalid format of basis set type in CPMD basis set file",1);
      /* number of shells */
      line = basis_cpmd_read_line(line,"number shells",file);
      if (sscanf(line,"%d%d",&n_shells,&n_prims)!=2)
        msg_error("invalid format of basis set type in CPMD basis set file",1);
      /* exponents & coefficients */
      line = basis_cpmd_read_line(line,"shell comment",file);
      line = basis_cpmd_read_line(line,"function exponents",file);
      de = str_parse_farray(line,&n);
      if (n_prims!=n)
        msg_error("inconsistent number of function exponents",1);
      id_p0 = id_p1 = 0;
      for (i=0; i<n_shells; i++) {
        line = basis_cpmd_read_line(line,"function coefficients",file);
        dc = str_parse_farray(line,&n);
        if (n_prims!=n)
          msg_error("inconsistent number of function coefficients",1);
        for (j=id_p0; j<n_prims; j++) {
          if (dc[j]==0.0)
            break;
          id_p1++;
          }
        s = basis_shell_new(1);
        s->type = basis_shell_id(type,1);
        s->n_prim = (id_p1-id_p0);
        s->exp = vec_falloc(s->n_prim);
        s->cf1 = vec_falloc(s->n_prim);
        for (k=0,j=id_p0; j<id_p1; j++,k++) {
          s->exp[k] = de[j];
          s->cf1[k] = dc[j];
          }
        list_add_end_p(st_s,s);
        vec_ffree(dc);
        id_p0 = id_p1;
        }
      vec_ffree(de);
      }
    /* convert list of shells to array */
    c->n_shells = st_s->num;
    c->shell = basis_shell_new(c->n_shells);
    for (id_s=0,p=st_s->first; p; p=p->l_next,id_s++) {
      s = (struct basis_shell*)p->l_data;
      c->shell[id_s].type = s->type;
      c->shell[id_s].n_prim = s->n_prim;
      c->shell[id_s].exp = s->exp;
      c->shell[id_s].cf1 = s->cf1;
      c->shell[id_s].cf2 = s->cf2;
      }
    list_free(st_s,NULL);
    /* store data in list */
    list_add_end_p(st_c,c);
    }
  /* close file */
  if (f && f[0])
    file_close(file);
  /* convert list of centers to array */
  b->n_centers = st_c->num;
  b->center = basis_center_new(b->n_centers);
  for (id_c=0,p=st_c->first; p; p=p->l_next,id_c++) {
    c = (struct basis_center*)p->l_data;
    b->center[id_c].type = c->type;
    b->center[id_c].n_shells = c->n_shells;
    b->center[id_c].shell = c->shell;
    }
  list_free(st_c,NULL);
  /* set global parameters */
  basis_def_parm(b);
  }

/* -------------------------------------------------------------------------- */

/* write basis set in CPMD format to file 
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_write_cpmd(struct basis *b, char *name) {
  unsigned ia,ic,id,ip,is,ir,in,n_shell,n_prim;
  unsigned at[120],sh[BASIS_SHELL_NUM],*cf_id;
  FILE *file = stdout;
  struct basis_center *c;
  struct basis_shell *s;
  if (!b)
    return;
  if (name && name[0])
    file = file_open(name,"w");
  /* select different atoms */
  vec_uset(at,0,120);
  for (ic=0; ic<b->n_centers; ic++)
    at[b->center[ic].type] = 1;
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      /* print atom name */
      fprintf(file,"%s\n",atom_name_full(ia));
      /* set shell center */
      ic = 0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      /* count shells */
      vec_uset(sh,0,BASIS_SHELL_NUM);
      for (is=0; is<c->n_shells; is++)
        sh[c->shell[is].type]++;
      if (sh[BASIS_SHELL_SP]) {
        sh[BASIS_SHELL_S]+=sh[BASIS_SHELL_SP];
        sh[BASIS_SHELL_P]+=sh[BASIS_SHELL_SP];
        sh[BASIS_SHELL_SP] = 0;
        }
      n_shell = 0;
      for (is=BASIS_SHELL_S; is<BASIS_SHELL_NUM; is++)
        if (sh[is] && is!=BASIS_SHELL_SP)
          n_shell++;
      fprintf(file,"%5d\n",n_shell);
      /* loop over shells */
      for (id=BASIS_SHELL_S; id<BASIS_SHELL_NUM; id++)
        if (sh[id]) {
          fprintf(file,"Shell: %s\n",basis_shell_name(id));
          /* count exponents */
          n_shell = n_prim = 0;
          for (is=0; is<c->n_shells; is++) {
            s = c->shell+is;
            if (s->type==id || (s->type==BASIS_SHELL_SP &&
               (id==BASIS_SHELL_S || id==BASIS_SHELL_P))) {
              n_prim += s->n_prim;
              n_shell++;
              }
            }
          fprintf(file,"%5d%5d\ncomment\n",n_shell,n_prim);
          cf_id = vec_ualloc(n_shell);
          /* exponents */
          n_shell = n_prim = 0;
          for (is=0; is<c->n_shells; is++) {
            s = c->shell+is;
            if (s->type==id || (s->type==BASIS_SHELL_SP &&
               (id==BASIS_SHELL_S || id==BASIS_SHELL_P))) {
              for (ip=0; ip<s->n_prim; ip++)
                fprintf(file,"%18.10e",s->exp[ip]);
              cf_id[n_shell] = n_prim;
              n_prim += s->n_prim;
              n_shell++;
              }
            }
          fprintf(file,"\n");
          /* coefficients */
          ir = 0;
          for (is=0; is<c->n_shells; is++) {
            s = c->shell+is;
            if (s->type==id || (s->type==BASIS_SHELL_SP &&
               (id==BASIS_SHELL_S || id==BASIS_SHELL_P))) {
              for (in=0; in<cf_id[ir]; in++)
                fprintf(file,"%18.10e",0.0);
              for (ip=0; ip<s->n_prim; ip++)
                fprintf(file,"%18.10e",
                  (s->type==BASIS_SHELL_SP && id==BASIS_SHELL_P ?
                   s->cf2[ip] : s->cf1[ip]));
              for (in=cf_id[ir]+s->n_prim; in<n_prim; in++)
                fprintf(file,"%18.10e",0.0);
              fprintf(file,"\n");
              ir++;
              }
            }
          /* clean memory */
          vec_ufree(cf_id);
          }
      }
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
