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
#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* split string to words according to given deliminer

   str - the string which is split
   del - the deliminer
   num - number of words read from the string */
char **str_split(char *str, char del, unsigned *num) {
  unsigned i,id,n;
  char *word,**v = NULL;
  struct queue *q;
  /* initialization */
  n = str_length(str);
  word = str_new(n);
  q = queue_alloc();
  /* split string to queue */
  for (id=0,n=str_length(str); id<n;) {
    id = str_read_word_s(str,word,del,id);
    queue_sadd(q,word);
    }
  /* convert queue to vector of strings */
  (*num) = q->num;
  v = (char**)malloc(q->num*sizeof(char*));
  if (!v) 
    msg_error("cannot allocate memory for spit-string array",1);
  for (i=0; i<(*num); i++) 
    v[i] = queue_get(q);
  /* clean memory */
  str_free(word);
  queue_free(q);
  return(v);
  }

/* split string to words according to given string deliminer

   s - the string which is split
   d - the string deliminer
   n - number of words read from the string */
char** str_split_s(char *s, char *d, unsigned *n) {
  char *word,**v = NULL;
  unsigned id,is,iw,ns,nd;
  short del;
  struct queue *q;
  /* initialization */
  is = 0;
  iw = 0;
  ns = str_length(s);
  nd = str_length(d);
  word = str_new(ns);
  q = queue_alloc();
  /* go through the string 's' */
  while (is<ns && s[is]!='\0') {
    /* check deliminer 'd' */
    for (id=0,del=0; id<nd; id++)
      if (s[is]==d[id]) {
        del = 1;
        break;
        }
    /* skip character */
    if (del) {
      if (iw) {
        word[iw] = '\0';
        queue_sadd(q,word);
        iw = 0;
        }
      is++;
      }
    /* save character */
    else {
      word[iw] = s[is];
      iw++;
      is++;
      }
    }
  /* add last word */
  if (iw) {
    word[iw] = '\0';
    queue_sadd(q,word);
    }
  /* convert queue to vector of strings */
  (*n) = q->num;
  v = (char**)malloc(q->num*sizeof(char*));
  if (!v) 
    msg_error("cannot allocate memory for spit-string array",1);
  for (iw=0; iw<(*n); iw++) 
    v[iw] = queue_get(q);
  /* clean memory */
  queue_free(q);
  return(v);
  }

/* split string to words according to given word deliminer

   s - the string which is split
   d - the string word deliminer
   n - number of words read from the string */
char** str_split_w(char *s, char *d, unsigned *n) {
  char *word,**v = NULL;
  unsigned id,ir,is,iw,ns,nd;
  short del;
  struct queue *q;
  /* initialization */
  is = 0;
  iw = 0;
  ns = str_length(s);
  nd = str_length(d);
  word = str_new(ns);
  q = queue_alloc();
  /* go through the string 's' */
  while (is<ns && s[is]!='\0') {
    /* check deliminer 'd' */
    for (id=0,ir=is,del=1; id<nd; id++,ir++)
      if (s[ir]!=d[id])
        del = 0;
    /* skip deliminer */
    if (del) {
      if (iw) {
        word[iw] = '\0';
        queue_sadd(q,word);
        iw = 0;
        }
      is = ir;
      del = 0;
      }
    /* save word */
    else {
      word[iw] = s[is];
      iw++;
      is++;
      }
    }
  /* add last word */
  if (iw) {
    word[iw] = '\0';
    queue_sadd(q,word);
    }
  /* convert queue to vector of strings */
  (*n) = q->num;
  v = (char**)malloc(q->num*sizeof(char*));
  if (!v) 
    msg_error("cannot allocate memory for spit-string array",1);
  for (iw=0; iw<(*n); iw++) 
    v[iw] = queue_get(q);
  /* clean memory */
  queue_free(q);
  return(v);
  }

/* -------------------------------------------------------------------------- */
