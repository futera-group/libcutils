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

#ifndef ZF_LIB_CMN_QUEUE_H
#define ZF_LIB_CMN_QUEUE_H

/* -------------------------------------------------------------------------- */

/* structs and types */
struct qdata {
  struct qdata *q_next;
  void *q_data;
  };

struct queue {
  long unsigned num;
  struct qdata *data;
  };

/* -------------------------------------------------------------------------- */

/* allocate new queue structure */
struct queue *queue_alloc(void);
/* destroy allocated queue structure */
void queue_free(struct queue*);
/* clean queue structure */
void queue_clean(struct queue*);

/* add new unspecified data to the queue */
void queue_add(struct queue*, void*);
/* add integer to the queue */
long unsigned queue_iadd(struct queue*, int);
/* add unsigned int to the queue */
long unsigned queue_uadd(struct queue*, unsigned);
/* add short integer to the queue */
long unsigned queue_siadd(struct queue*, short);
/* add short unsigned integer to the queue */
long unsigned queue_suadd(struct queue*, short unsigned);
/* add long integer to the queue */
long unsigned queue_liadd(struct queue*, long);
/* add long unsigned integer to the queue */
long unsigned queue_luadd(struct queue*, long unsigned);
/* add double precision real to the queue */
long unsigned queue_fadd(struct queue*, double);
/* add long double precision real to the queue */
long unsigned queue_lfadd(struct queue*, long double);
/* add character to the queue */
long unsigned queue_cadd(struct queue*, char);
/* add string to the queue */
long unsigned queue_sadd(struct queue*, char*);

/* return pointer to the first data item in the queue  */
void *queue_get(struct queue*);
/* get integer from the queue */
long unsigned queue_iget(struct queue*, int*);
/* get unsigned int from the queue */
long unsigned queue_uget(struct queue*, unsigned*);
/* get short integer from the queue */
long unsigned queue_siget(struct queue*, short*);
/* get short unsigned integer from the queue */
long unsigned queue_suget(struct queue*, short unsigned*);
/* get long integer from the queue */
long unsigned queue_liget(struct queue*, long*);
/* get long unsigned integer from the queue */
long unsigned queue_luget(struct queue*, long unsigned*);
/* get double precision real from the queue */
long unsigned queue_fget(struct queue*, double*);
/* get long double precision real from the queue */
long unsigned queue_lfget(struct queue*, long double*);
/* get integer from the queue */
long unsigned queue_cget(struct queue*, char*);
/* get string from the queue */
long unsigned queue_sget(struct queue*, char*);

/* convert data saved in a queue to vector */
void *queue_conv_vec(struct queue*, short, long unsigned*);

/* -------------------------------------------------------------------------- */

#endif
