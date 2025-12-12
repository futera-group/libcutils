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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/print.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "qmc/basis.h"
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

#define NBO_SHELL_S  1
#define NBO_SHELL_P  2
#define NBO_SHELL_D  3
#define NBO_SHELL_F  4

/* -------------------------------------------------------------------------- */

/* read block of values from NBO basis set specification
 
   n    - number of expected values
   file - pointer to open input file stream */
double* basis_nbo_read_array(unsigned n, FILE *file) {
  char *line;
  unsigned i,j,id,nr,nv;
  double *v;
  id = 0;
  v = vec_falloc(n);
  nr = (n%4 ? n/4+1 : n/4);
  for (i=0; i<nr; i++) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of input file while reading data array",1);
    nv = (i+1<nr ? 4 : (n%4 ? n%4 : 4));
    for (j=0; j<nv; j++)
      if (sscanf(line+2+18*j,"%lf",&(v[id++]))!=1)
        msg_error("invalid format of data array",1);
    str_free(line);
    }
  return(v);
  }

/* read basis set in NBO format from file
 
   b - pointer to basis set data struct
   f - name of the file */
void basis_nbo_read(struct basis *b, char *f) {
  unsigned i,j,k,n,fn,*fi,**id,a[3];
  char *line;
  double *v,norm;
  FILE *file = NULL;
  if (f && f[0])
    file = file_open(f,"r");
  /* header */
  for (i=0; i<3; i++) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of \"%s\" file while reading header",1,f);
    str_free(line);
    }
  /* number of centers and shells */
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading numbers",1,f);
  if (sscanf(line,"%d%d%d",&(b->n_centers),&(b->n_cont_shells),
      &(b->n_prim_shells))!=3)
    msg_error_f("invalid format of number of centers and shells",1);
  line = str_free(line);
  b->center=basis_center_new(b->n_centers);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading numbers",1,f);
  line = str_free(line);
  /* centers */
  for (i=0; i<b->n_centers; i++) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of \"%s\" file while reading centers",1,f);
    if (sscanf(line,"%d%lf%lf%lf",&(b->center[i].type),
       &(b->center[i].coord[0]),&(b->center[i].coord[1]),
       &(b->center[i].coord[2]))!=4)
      msg_error_f("invalid format of center type and coordinates",1);
    line = str_free(line);
    }
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading centers",1,f);
  line = str_free(line);
  /* shell types */
  id = mat_ualloc(b->n_cont_shells,100);
  for (i=0; i<b->n_cont_shells; i++) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of \"%s\" file while reading shell num",1,f);
    if (sscanf(line,"%d%d%d%d",&(id[i][0]),&(id[i][1]),&(id[i][2]),
        &(id[i][3]))!=4)
      msg_error_f("invalid format of shell numbers",1);
    line = str_free(line);
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of \"%s\" file while reading shell ID",1,f);
    fi = str_parse_uarray(line,&fn);
    if (id[i][1]!=fn)
      msg_error_f("invalid number of shell IDs",1);
    for (j=0; j<fn; j++)
      id[i][4+j]=fi[j];
    vec_ufree(fi);
    line = str_free(line);
    }
  /* convert ID array to basis set struct */
  for (i=0; i<b->n_centers; i++) {
    /* number of shells */
    n = 0;
    for (j=0; j<b->n_cont_shells; j++)
      if (id[j][0]==(i+1))
        n++;
    b->center[i].n_shells = n;
    b->center[i].shell = basis_shell_new(n);
    /* shell type */
    n = 0;
    for (j=0; j<b->n_cont_shells; j++)
      if (id[j][0]==(i+1)) {
        switch (id[j][1]) {
          case  1: b->center[i].shell[n].type=BASIS_SHELL_S;  break;
          case  3: b->center[i].shell[n].type=BASIS_SHELL_P;  break;
          case  4: b->center[i].shell[n].type=BASIS_SHELL_SP; break;
          case  5: b->center[i].shell[n].type=BASIS_SHELL_Dp; break;
          case  6: b->center[i].shell[n].type=BASIS_SHELL_Dc; break;
          case  7: b->center[i].shell[n].type=BASIS_SHELL_Fp; break;
          case 10: b->center[i].shell[n].type=BASIS_SHELL_Fc; break;
          }
        b->center[i].shell[n].n_prim = id[j][3];
        b->center[i].shell[n].exp = vec_falloc(id[j][3]);
        b->center[i].shell[n].cf1 = vec_falloc(id[j][3]);
        if (b->center[i].shell[n].type==BASIS_SHELL_SP)
          b->center[i].shell[n].cf2 = vec_falloc(id[j][3]);
        n++;
        }
    }
  mat_ufree(id,b->n_cont_shells);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading shell types",1,f);
  line = str_free(line);
  /* exponents */
  v = basis_nbo_read_array(b->n_prim_shells,file);
  n = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++)
         b->center[i].shell[j].exp[k] = v[n++];
  vec_ffree(v);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading exponents",1,f);
  line = str_free(line);
  /* coefficients */
  v = basis_nbo_read_array(b->n_prim_shells,file);
  n = 0;
  a[0] = a[1] = a[2] = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
         if (b->center[i].shell[j].type==BASIS_SHELL_S ||
             b->center[i].shell[j].type==BASIS_SHELL_SP) {
           norm = gto_norm(b->center[i].shell[j].exp[k],a);
           b->center[i].shell[j].cf1[k] = v[n]/norm;
           }
         n++;
         }
  vec_ffree(v);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading s-coeff",1,f);
  line = str_free(line);
  v = basis_nbo_read_array(b->n_prim_shells,file);
  n = 0;
  a[0] = 1;
  a[1] = a[2] = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
         if (b->center[i].shell[j].type==BASIS_SHELL_P) {
           norm = gto_norm(b->center[i].shell[j].exp[k],a);
           b->center[i].shell[j].cf1[k] = v[n]/norm;
           }
         else if (b->center[i].shell[j].type==BASIS_SHELL_SP) {
           norm = gto_norm(b->center[i].shell[j].exp[k],a);
           b->center[i].shell[j].cf2[k] = v[n]/norm;
           }
         n++;
         }
  vec_ffree(v);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading p-coeff",1,f);
  line = str_free(line);
  v = basis_nbo_read_array(b->n_prim_shells,file);
  n = 0;
  a[0] = 2;
  a[1] = a[2] = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
         if (b->center[i].shell[j].type==BASIS_SHELL_Dp ||
             b->center[i].shell[j].type==BASIS_SHELL_Dc) {
           norm = gto_norm(b->center[i].shell[j].exp[k],a);
           b->center[i].shell[j].cf1[k] = v[n]/norm;
           }
         n++;
         }
  vec_ffree(v);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of \"%s\" file while reading d-coeff",1,f);
  line = str_free(line);
  v = basis_nbo_read_array(b->n_prim_shells,file);
  n = 0;
  a[0] = 3;
  a[1] = a[2] = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
         if (b->center[i].shell[j].type==BASIS_SHELL_Fp ||
             b->center[i].shell[j].type==BASIS_SHELL_Fc) {
           norm = gto_norm(b->center[i].shell[j].exp[k],a);
           b->center[i].shell[j].cf1[k] = v[n]/norm;
           }
         n++;
         }
  vec_ffree(v);
  /* close file */
  if (f && f[0])
    file_close(file);
  /* set global parameters */
  basis_def_parm(b);
  }

/* -------------------------------------------------------------------------- */

/* write basis shell info in NBO format
 
   b    - pointer to basis set struct
   file - pointer to open output stream */
void basis_nbo_write_shell(struct basis *b, FILE *file) {
  unsigned i,j,k,n = 0,x = 0,n_prim = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++) {
      fprintf(file,"%7d%6d%6d%6d\n",i+1,b->center[i].shell[j].n_bfce,
        n_prim+1,b->center[i].shell[j].n_prim);
      switch (b->center[i].shell[j].type) {
        case BASIS_SHELL_S:  n =  1; x =   0; break;
        case BASIS_SHELL_P:
        case BASIS_SHELL_SP: n =  3; x = 100; break;
        case BASIS_SHELL_Dp: n =  5; x = 250; break;
        case BASIS_SHELL_Dc: n =  6; x = 200; break;
        case BASIS_SHELL_Fp: n =  7; x = 350; break;
        case BASIS_SHELL_Fc: n = 10; x = 300; break;
        default:
          fprintf(file,"\n");
          break;
        }
      fprintf(file," ");
      if (b->center[i].shell[j].type==BASIS_SHELL_SP)
        fprintf(file,"%6d",1);
      for (k=1; k<=n; k++)
        fprintf(file,"%6d",x+k);
      fprintf(file,"\n");
      n_prim += b->center[i].shell[j].n_prim;
      }
  }

/* write basis fce exponents in NBO format
 
   b    - pointer to basis set struct
   file - pointer to open output stream */
void basis_nbo_write_exp(struct basis *b, FILE *file) {
  unsigned i,j,k,n;
  n = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
        if (!n || (n%4)==0)
          fprintf(file,"  ");
        fprintf(file,"%18.9E",b->center[i].shell[j].exp[k]);
        if ((++n%4)==0)
          fprintf(file,"\n");
        }
  if (n%4)
    fprintf(file,"\n");
  }

/* write basis fce coefficients in NBO format
 
   b    - pointer to basis set struct
   file - pointer to open output stream
   t    - type of the shell (s,p,d,f) */
void basis_nbo_write_coeff(struct basis *b, FILE *file, short t) {
  unsigned i,j,k,n,a[3];
  double coeff,norm;
  short calc;
  n = 0;
  for (i=0; i<b->n_centers; i++)
    for (j=0; j<b->center[i].n_shells; j++)
      for (k=0; k<b->center[i].shell[j].n_prim; k++) {
        calc = 0;
        a[0] = a[1] = a[2] = 0;
        switch (t) {
          case NBO_SHELL_S:
            if (b->center[i].shell[j].type==BASIS_SHELL_S ||
                b->center[i].shell[j].type==BASIS_SHELL_SP)
              calc = 1;
            break;
          case NBO_SHELL_P:
            if (b->center[i].shell[j].type==BASIS_SHELL_P ||
                b->center[i].shell[j].type==BASIS_SHELL_SP) {
              a[0] = 1;
              calc = 1;
              }
            break;
          case NBO_SHELL_D:
            if (b->center[i].shell[j].type==BASIS_SHELL_Dc) {
              a[0] = 2;
              calc = 1;
              }
            break;
          case NBO_SHELL_F:
            if (b->center[i].shell[j].type==BASIS_SHELL_Fc) {
              a[0] = 3;
              calc = 1;
              }
            break;
          }
        if (calc) {
          norm = gto_norm(b->center[i].shell[j].exp[k],a);
          coeff = norm*((t==NBO_SHELL_P && 
            b->center[i].shell[j].type==BASIS_SHELL_SP) ? 
            b->center[i].shell[j].cf2[k] :
            b->center[i].shell[j].cf1[k]);
          }
        else
          coeff = 0.0;
        if (!n || (n%4)==0)
          fprintf(file,"  ");
        fprintf(file,"%18.9E",coeff);
        if ((++n%4)==0)
          fprintf(file,"\n");
        }
  if (n%4)
    fprintf(file,"\n");
  }

/* write basis set in NBO format to file 
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_nbo_write(struct basis *b, char *name) {
  unsigned i,j;
  FILE *file = stdout;
  if (!b)
    return;
  if (name && name[0])
    file = file_open(name,"w");
  fprintf(file," Title\n");
  fprintf(file," Basis set information needed for plotting orbitals\n ");
  print_fHline(file,'-',75);
  fprintf(file,"%7d%6d%6d\n ",b->n_centers,b->n_cont_shells,b->n_prim_shells);
  print_fHline(file,'-',75);
  /* centers */
  for (i=0; i<b->n_centers; i++) {
    fprintf(file,"%5d",b->center[i].type);
    for (j=0; j<3; j++)
      fprintf(file,"%14.9f",b->center[i].coord[j]);
    fprintf(file,"\n");
    }
  fprintf(file," ");
  print_fHline(file,'-',75);
  /* shell info */
  basis_nbo_write_shell(b,file);
  fprintf(file," ");
  print_fHline(file,'-',75);
  /* exponents */
  basis_nbo_write_exp(b,file);
  /* coefficients  */
  fprintf(file,"\n");
  basis_nbo_write_coeff(b,file,NBO_SHELL_S);
  fprintf(file,"\n");
  basis_nbo_write_coeff(b,file,NBO_SHELL_P);
  fprintf(file,"\n");
  basis_nbo_write_coeff(b,file,NBO_SHELL_D);
  fprintf(file,"\n");
  basis_nbo_write_coeff(b,file,NBO_SHELL_F);
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
