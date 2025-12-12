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

#include <stdlib.h>
#include <cmn/file.h>
#include <cmn/list.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* read one non empty line from basis set file in gaussian format
 
   line - storage for the line
   file - pointer to open file stream */
char* basis_gauss_read_line(char *line, FILE *file) {
  for (line=str_free(line),line=str_read_line_new(file); line;
       line=str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (str_length(line)<2)
      return(NULL);
    if (str_length(line)>2 && line[0]!='!') {
      return(line);
      }
    }
  return(line);
  }

/* read basis set in gaussian format from open file
 
   b    - pointer to basis set data struct
   file - open stream of the file */
void basis_gauss_read_f(struct basis *b, FILE *file) {
  char name[256],type[256],*line = NULL;
  unsigned ip,id_s,id_c;
  struct basis_center *c;
  struct basis_shell *s;
  struct list *st_c,*st_s;
  struct ldata *p;
  /* read data from file */
  st_c = list_alloc();
  for (line=basis_gauss_read_line(line,file); line;
       line=basis_gauss_read_line(line,file)) {
    /* new center */
    c = basis_center_new(1);
    if (sscanf(line,"%s",name)!=1)
      msg_error("invalid format of basis set header",1);
    c->type = atom_num(name);
    /* read shells */
    st_s = list_alloc();
    for (line=basis_gauss_read_line(line,file); line && line[0]!='*';
         line=basis_gauss_read_line(line,file)) {
      s = basis_shell_new(1);
      if (sscanf(line,"%s%d",type,&(s->n_prim))!=2)
        msg_error_f("invalid format of shell in %s basis set",1,name);
      s->type = basis_shell_id(type,0);
      /* read primitive shells */
      s->exp = vec_falloc(s->n_prim);
      s->cf1 = vec_falloc(s->n_prim);
      if (s->type==BASIS_SHELL_SP)
        s->cf2 = vec_falloc(s->n_prim);
      for (ip=0; ip<s->n_prim; ip++) {
        line = basis_gauss_read_line(line,file);
        if (!line)
          msg_error_f("invalid format of %s basis set, %s shell",1,name,type);
        if (s->type==BASIS_SHELL_SP) {
          if (sscanf(line,"%lf%lf%lf",&(s->exp[ip]),&(s->cf1[ip]),
            &(s->cf2[ip]))!=3)
            msg_error_f("invalid format of %s basis set, %s shell, prim %d",1,
              name,type,ip+1);
            }
        else {
          if (sscanf(line,"%lf%lf",&(s->exp[ip]),&(s->cf1[ip]))!=2)
            msg_error_f("invalid format of %s basis set, %s shell, prim %d",1,
              name,type,ip+1);
            }
        line = str_free(line);
        }
      list_add_end_p(st_s,s);
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

/* read basis set in gaussian format from file
 
   b - pointer to basis set data struct
   f - name of the file */
void basis_gauss_read(struct basis *b, char *f) {
  FILE *file = NULL;
  if (f && f[0])
    file = file_open(f,"r");
  basis_gauss_read_f(b,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */

/* write basis set in gaussian format to file 
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_gauss_write(struct basis *b, char *name) {
  unsigned i,j,k,m,at[120];
  FILE *file = stdout;
  struct basis_shell *bs;
  if (!b)
    return;
  if (name && name[0])
    file = file_open(name,"w");
  /* select different atoms */
  vec_uset(at,0,120);
  for (i=0; i<b->n_centers; i++)
    at[b->center[i].type] = 1;
  for (i=0; i<120; i++)
    if (at[i]) {
      /* print atom name */
      fprintf(file,"%s  0\n",atom_name(i));
      j = 0;
      while (b->center[j].type!=i)
        j++;
      for (k=0; k<b->center[j].n_shells; k++) {
        bs = b->center[j].shell+k;
        /* print shell type and number of primitives */
        fprintf(file,"%-2s %3d  1.0\n",basis_shell_name(bs->type),bs->n_prim);
        for (m=0; m<bs->n_prim; m++) {
          /* print coefficients and exponents of primitives */
          fprintf(file,"%18.10f %18.10f",bs->exp[m],bs->cf1[m]);
          if (bs->type==BASIS_SHELL_SP)
            fprintf(file," %18.10f",bs->cf2[m]);
          fprintf(file,"\n");
          }
        }
      fprintf(file,"****\n");
      }
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
