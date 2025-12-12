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
#include <string.h>
#include "cmn/config.h"
#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"
#include "cmn/types.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* assign data value from keyword line with specified type to variable

   s - keyword line from the config file
   k - name of the keyword
   t - type of the data
   d - pointer to data storage */
void cfg_assign_rank0(char *s, char *k, short t, void *d) {
  char key[1024],val[1024];
  if (sscanf(s,"%s%s",key,val)!=2)
    msg_error_f("invalid format of config line with keyword \"%s\"",1,k);
  if (!type_read(val,d,t))
    msg_error_f("cannot assign scalar value to keyword \"%s\"",1,k);
  }

/* assign data value from keyword line with specified type to variable

   q - storage of the data lines read from file 
   k - name of the keyword
   t - type of the data
   n - number of expected data elements
   d - pointer to data storage */
void cfg_assign_rank1(struct queue *q, char *k, short t, unsigned n, void *d) {
  char *line;
  unsigned i,id,*vi,nv,*ti;
  double *vf,*tf;
  long double *vx,*tx;
  switch (t) {
    case TYPE_UINT: /* unsigned integer */
      id = 0;
      ti = *((unsigned**)d);
      while (q->num) {
        line = queue_get(q);
        if (line) {
          if (strstr(line,",") || strstr(line,"-"))
            vi = str_parse_unum(line,&nv);
          else 
            vi = str_parse_uarray(line,&nv);
          if ((id+nv)>n)
            msg_error_f("number of rank-1 data values of \"%s\""
              " is %d, not %d as specified",1,k,nv,n);
          for (i=0; i<nv; i++)
            ti[id++] = vi[i];
          }
        str_free(line);
        } 
      break;
    case TYPE_DOUBLE: /* double */
      id = 0;
      tf = *((double**)d);
      while (q->num) {
        line = queue_get(q);
        if (line) {
          vf = str_parse_farray(line,&nv);
          if ((id+nv)>n)
            msg_error_f("number of rank-1 data values of \"%s\""
              " is %d, not %d as specified",1,k,nv,n);
          for (i=0; i<nv; i++)
            tf[id++] = vf[i];
          }
        str_free(line);
        }
      break;
    case TYPE_LDOUBLE: /* long double */
      id = 0;
      tx = *((long double**)d);
      while (q->num) {
        line = queue_get(q);
        if (line) {
          vx = str_parse_lfarray(line,&nv);
          if ((id+nv)>n)
            msg_error_f("number of rank-1 data values of \"%s\""
              " is %d, not %d as specified",1,k,nv,n);
          for (i=0; i<nv; i++)
            tx[id++] = vx[i];
          }
        str_free(line);
        }
      break;
    default:
      msg_error_f("unsupported rank-1 data type (%d) of \"%s\" keyword",1,t,k);
    } 
  } 

/* assign data value from keyword line with specified type to variable

   q - storage of the data lines read from file 
   k - name of the keyword
   t - type of the data
   n - number of expected data elements
   m - storage for rank-1 array dimensions
   d - pointer to data storage */
void cfg_assign_rank2(struct queue **q, char *k, short t, unsigned n,
  unsigned *m, void *d) {
  unsigned i,j,id,*vu,nv;
  double *vf;
  long double *vx;
  char *line;
  struct queue *x;
  switch (t) {
    case TYPE_UINT: /* unsigned integer */
      x = queue_alloc();
      /* rank-1 arrays */
      for (i=0; i<n; i++) {
        /* parse saved data lines */
        while (q[i]->num) {
          line = queue_get(q[i]);
          if (line) {
            if (strstr(line,",") || strstr(line,"-"))
              vu = str_parse_unum(line,&nv);
            else 
              vu = str_parse_uarray(line,&nv);
            /* save numbers */
            for (j=0; j<nv; j++)
              queue_uadd(x,vu[j]);
            vu = vec_ufree(vu);
            }
          str_free(line);
          } 
        /* save number to array */
        id = 0;
        m[i] = x->num;
        (*(unsigned***)d)[i] = vec_ualloc(m[i]);
        while (x->num)
          queue_uget(x,(*(unsigned***)d)[i]+id++);
        }
      /* clean memory */
      queue_free(x);
      break;
    case TYPE_DOUBLE: /* double-precision real */
      x = queue_alloc();
      /* rank-1 arrays */
      for (i=0; i<n; i++) {
        /* parse saved data lines */
        while (q[i]->num) {
          line = queue_get(q[i]);
          if (line) {
            vf = str_parse_farray(line,&nv);
            for (j=0; j<nv; j++)
              queue_fadd(x,vf[j]);
            vf = vec_ffree(vf);
            }
          str_free(line);
          } 
        /* save number to array */
        id = 0;
        m[i] = x->num;
        (*(double***)d)[i] = vec_falloc(m[i]);
        while (x->num)
          queue_fget(x,(*(double***)d)[i]+id++);
        }
      /* clean memory */
      queue_free(x);
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      x = queue_alloc();
      /* rank-1 arrays */
      for (i=0; i<n; i++) {
        /* parse saved data lines */
        while (q[i]->num) {
          line = queue_get(q[i]);
          if (line) {
            vx = str_parse_lfarray(line,&nv);
            for (j=0; j<nv; j++)
              queue_lfadd(x,vx[j]);
            vx = vec_lffree(vx);
            }
          str_free(line);
          } 
        /* save number to array */
        id = 0;
        m[i] = x->num;
        (*(long double***)d)[i] = vec_lfalloc(m[i]);
        while (x->num)
          queue_lfget(x,(*(long double***)d)[i]+id++);
        }
      /* clean memory */
      queue_free(x);
      break;
    default:
      msg_error_f("unsupported rank-2 data type (%d) of \"%s\" keyword",1,t,k);
    }
  }

/* -------------------------------------------------------------------------- */
