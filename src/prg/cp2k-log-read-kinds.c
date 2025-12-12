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
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* save primitive funcions read from the log file to orbital data struct
 
   qn - atomic orbital quantum numbers
   qa - primitive function exponents
   qc - primitive function coefficients
   qs - shell where the AO is to be added */
void cp2k_log_read_ao(int *qn, struct queue *qa, struct queue *qc,
  struct queue *qs) {
  unsigned i;
  struct cp2k_ao *ao;
  /* sanity check */
  if (qa->num!=qc->num)
    msg_error_f("inconsistent number of primitive function exponents (%ld)"
      " and coefficients (%ld)",1,qa->num,qc->num);
  if (qa->num>0) {
    /* memory allocation */
    ao = cp2k_ao_new(1);
    /* quantum numbers */
    ao->qnum[0] = qn[0];
    ao->qnum[1] = qn[1];
    ao->qnum[2] = qn[2];
    /* primitive functions */
    ao->n_prim_fce = qa->num;
    ao->exp = vec_falloc(ao->n_prim_fce);
    ao->coeff = vec_falloc(ao->n_prim_fce);
    for (i=0; i<ao->n_prim_fce; i++) {
      queue_fget(qa,&(ao->exp[i]));
      queue_fget(qc,&(ao->coeff[i]));
      }
    /* save AO to shell */
    queue_add(qs,ao);
    }
  }

/* read atomic kinds from cp2k output file
 
   d - pointer to cp2k data struct
   f - open file stream */
void cp2k_log_read_kinds(struct cp2k_dat *d, FILE *f) {
  char *line,flag[80],**s;
  unsigned i,j,k,n,id,ns,nc;
  int m,tp[3];
  double a,c;
  struct cp2k_ao *ao;
  struct queue *qa,*qc,*qs;
  /* find the atom data block */
  rewind(f);
  line = str_ffind_b_new(f," ATOMIC KIND INFORMATION");
  if (line) {
    line = str_free(line);
    /* storage memory */
    qa = queue_alloc();
    qc = queue_alloc();
    qs = queue_alloc();
    /* read kind data one-by-one */
    if (!d->kind)
      d->kind = cp2k_kind_new(d->n_atom_kinds);
    for (i=0; i<d->n_atom_kinds; i++) {
      sprintf(flag,"%3d. Atomic kind:",i+1);
      line = str_ffind_b_new(f,flag);
      if (!line)
        msg_error_f("cp2k atom kind #%d block not found",1,i+1);
      /* number of atoms */
      str_trim(line);
      for (n=str_length(line); n>0 && line[n-1]!=' '; n--);
      if (sscanf(line+n,"%d",&(d->kind[i].n_atoms))!=1)
        msg_error_f("cannot read number of atoms in atom kind #%d block",1,i+1);
      line = str_free(line);
      /* number of orbitals */
      for (j=0; j<9; j++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of cp2k log file",1);
        if (j>2) {
          if (sscanf(line+70,"%d",&m)!=1)
            msg_error_f("cannot read orbital numbers of atom kind #%d",1,i+1);
          switch (j) {
            case 3: d->kind[i].n_shell_sets = m; break;
            case 4: d->kind[i].n_shells = m;     break;
            case 5: d->kind[i].n_prim_fce = m;   break;
            case 6: d->kind[i].n_cart_fce = m;   break;
            case 7: d->kind[i].n_sphr_fce = m;   break;
            case 8: d->kind[i].norm_type = m;    break;
            }
          }
        line = str_free(line);
        }
      /* atomic orbitals */
      for (j=0; j<5; j++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of cp2k log file",1);
        line = str_free(line);
        }
      if (!d->kind[i].shell)
        d->kind[i].shell = cp2k_shell_new(d->kind[i].n_shells);
      for (j=0; j<d->kind[i].n_shells; j++) {
        id = tp[0] = tp[1] = tp[2] = 0;
        /* shell block */
        for (;;) {
          line = str_read_line_new(f);
          if (!line)
            msg_error("unexpected end of cp2k log file",1);
          if (str_length(line)<10)
            break;
          s = str_split(line,' ',&n);
          line = str_free(line);
          switch (n) {
            case 2: /* the same sub-type */
              if (sscanf(s[0],"%lf",&a)!=1)
                msg_error_f("cannot read orbital exponent of"
                  " atom kind #%d, shell #%d",1,i+1,j+1);
              if (sscanf(s[1],"%lf",&c)!=1)
                msg_error_f("cannot read orbital coefficient of"
                  " atom kind #%d, shell #%d",1,i+1,j+1);
              break;
            case 5: /* new orbital sub-type */
              /* save AO data */
              cp2k_log_read_ao(tp,qa,qc,qs);
              /* read next AO */
              if (sscanf(s[0],"%d",&id)!=1)
                msg_error_f("cannot read shell-set ID of"
                  " atom kind #%d, shell #%d",1,i+1,j+1);
              if (id<1 || id>d->kind[i].n_shell_sets)
                msg_error_f("invalid shell-set ID (%d) of"
                  " atom kind #%d, shell #%d",1,id,i+1,j+1);
              cp2k_orbital_type_id(s[2],&(tp[0]),&(tp[1]),&(tp[2]));
              if (sscanf(s[3],"%lf",&a)!=1)
                msg_error_f("cannot read orbital exponent of"
                  " atom kind #%d, shell #%d",1,i+1,j+1);
              if (sscanf(s[4],"%lf",&c)!=1)
                msg_error_f("cannot read orbital coefficient of"
                  " atom kind #%d, shell #%d",1,i+1,j+1);
              break;
            default:
              msg_error_f("unknown format of atom kind #%d"
                " orbital block",1,i+1);
            }
          vec_sfree(s,n);
          /* save primitive functions */
          queue_fadd(qa,a);
          queue_fadd(qc,c);
          }
        /* save AO data */
        cp2k_log_read_ao(tp,qa,qc,qs);
        /* save shell data */
        d->kind[i].shell[j].set = id;
        d->kind[i].shell[j].type = 0;
        d->kind[i].shell[j].n_cart_fce = qs->num;
        d->kind[i].shell[j].ao = cp2k_ao_new(d->kind[i].shell[j].n_cart_fce);
        for (k=0; k<d->kind[i].shell[j].n_cart_fce; k++) {
          ao = queue_get(qs);
          vec_icopy(d->kind[i].shell[j].ao[k].qnum,ao->qnum,3);
          d->kind[i].shell[j].ao[k].n_prim_fce = 
            ao->n_prim_fce;
          d->kind[i].shell[j].ao[k].exp = 
            vec_fcopy_new(ao->exp,ao->n_prim_fce);
          d->kind[i].shell[j].ao[k].coeff = 
            vec_fcopy_new(ao->coeff,ao->n_prim_fce);
          if (k==0)
            d->kind[i].shell[j].type = ao->qnum[1];
          else if (d->kind[i].shell[j].type!=ao->qnum[1])
            msg_error_f("more than one basis-function types in shell #%d of "
              "atom kind #%d",1,i+1,j+1);
          }
        d->kind[i].shell[j].n_sphr_fce = 2*d->kind[i].shell[j].type+1;
        }
      }
    /* check consistency */
//    ns = nc = 0;
//    for (i=0; i<d->n_atoms; i++) {
//      nc += d->kind[d->atom[i].kind].n_cart_fce;
//      ns += d->kind[d->atom[i].kind].n_sphr_fce;
//      }
//    if (d->n_sphr_fce!=ns)
//      msg_error_f("inconsistent number of spherical functions (%d/%d)"
//        " in basis set",1,d->n_sphr_fce,ns);
//    if (d->n_cart_fce!=nc)
//      msg_error_f("inconsistent number of cartesian functions (%d/%d)"
//        " in basis set",1,d->n_cart_fce,nc);
    /* clean memory */
    queue_free(qa);
    queue_free(qc);
    queue_free(qs);
    }
  }

/* -------------------------------------------------------------------------- */
