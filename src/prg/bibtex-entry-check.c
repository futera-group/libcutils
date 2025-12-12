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
#include <cmn/queue.h>
#include <cmn/string.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* check mandatory options and print out missing ones 
 
   e - pointer to bibtex entry data struct */
void bibtex_entry_check(struct bib_entry *e) {
  char *s;
  unsigned n = 0;
  struct queue *qm,*qo;
  /* mandatory and optional fields */
  qm = queue_alloc();
  qo = queue_alloc();
  /* entry types */
  switch (e->type) {
    case BIB_TYPE_ARTICLE: /* article */
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->journal || !e->journal[0])
        queue_sadd(qm,"journal");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->volume)
        queue_sadd(qm,"volume");
      if (!e->number)
        queue_sadd(qo,"number");
      if (!e->page[0])
        queue_sadd(qo,"pages");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_BOOK: /* book */
      if ((!e->author || !e->n_authors) && (!e->editor || !e->n_editors))
        queue_sadd(qm,"author/editor");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->publisher || !e->publisher[0])
        queue_sadd(qm,"publisher");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->number)
        queue_sadd(qo,"number");
      if (!e->page[0])
        queue_sadd(qo,"pages");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_BOOKLET: /* booklet */
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->author || !e->n_authors)
        queue_sadd(qo,"author");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->year)
        queue_sadd(qo,"year");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_CONF: /* conference proceeding */
    case BIB_TYPE_INPROC:
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->book_title || !e->book_title[0])
        queue_sadd(qm,"book_title");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->editor || !e->n_editors)
        queue_sadd(qo,"editor");
      if (!e->volume && !e->number)
        queue_sadd(qo,"volume/number");
      if (!e->page[0])
        queue_sadd(qo,"pages");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->organization || !e->organization[0])
        queue_sadd(qo,"organization");
      if (!e->publisher || !e->publisher[0])
        queue_sadd(qo,"publisher");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_INBOOK: /* part of book */
      if ((!e->author || !e->n_authors) && (!e->editor || !e->n_editors))
        queue_sadd(qm,"author/editor");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->chapter || !e->chapter[0])
        queue_sadd(qm,"chapter");
      if (!e->page[0])
        queue_sadd(qm,"pages");
      if (!e->publisher || !e->publisher[0])
        queue_sadd(qm,"publisher");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->volume && !e->number)
        queue_sadd(qo,"volume/number");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_INCOLL: /* titled part of book */
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->book_title || !e->book_title[0])
        queue_sadd(qm,"book_title");
      if (!e->publisher || !e->publisher[0])
        queue_sadd(qm,"publisher");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->editor || !e->n_editors)
        queue_sadd(qo,"editor");
      if (!e->volume && !e->number)
        queue_sadd(qo,"volume/number");
      if (!e->chapter || !e->chapter[0])
        queue_sadd(qo,"chapter");
      if (!e->page[0])
        queue_sadd(qo,"pages");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_MANUAL: /* technical documentation */
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->author || !e->n_authors)
        queue_sadd(qo,"author");
      if (!e->organization || !e->organization[0])
        queue_sadd(qo,"organization");
      if (!e->address[0])
        queue_sadd(qo,"address");
      if (!e->year)
        queue_sadd(qo,"year");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_MTHESIS: /* master thesis */
    case BIB_TYPE_PHDTHESIS: /* ph.d. thesis */
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->school || !e->school[0])
        queue_sadd(qm,"school");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_MISC: /* others */
      if (!e->author || !e->n_authors)
        queue_sadd(qo,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qo,"title");
      if (!e->year)
        queue_sadd(qo,"year");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_PROC: /* proceedings */
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->editor || !e->n_editors)
        queue_sadd(qo,"editor");
      if (!e->volume && !e->number)
        queue_sadd(qo,"volume/number");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->publisher || !e->publisher[0])
        queue_sadd(qo,"publisher");
      if (!e->organization || !e->organization[0])
        queue_sadd(qo,"organization");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_TECHREP: /* technical report */
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->organization || !e->organization[0])
        queue_sadd(qm,"organization");
      if (!e->year)
        queue_sadd(qm,"year");
      if (!e->number)
        queue_sadd(qo,"number");
      if (!e->address || !e->address[0])
        queue_sadd(qo,"address");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    case BIB_TYPE_UNPUB: /* unpublished */
      if (!e->author || !e->n_authors)
        queue_sadd(qm,"author");
      if (!e->title || !e->title[0])
        queue_sadd(qm,"title");
      if (!e->year)
        queue_sadd(qo,"year");
      if (!e->keyword || !e->n_keywords)
        queue_sadd(qo,"keywords");
      break;
    }
  /* missing mandatory fields */
  if (qm->num) {
    printf("\n");
    printf("Missing mandatory fields: ");
    n = 0;
    while (qm->num) {
      s = queue_get(qm);
      if (n)
        printf(", ");
      printf("%s",s);
      str_free(s);
      n++;
      }
    printf("\n");
    }
  /* missing optional fields */
  if (qo->num) {
    if (!n)
      printf("\n");
    printf("Missing optional fields:  ");
    n = 0;
    while (qo->num) {
      s = queue_get(qo);
      if (n)
        printf(", ");
      printf("%s",s);
      str_free(s);
      n++;
      }
    printf("\n");
    }
  /* clean memory */
  queue_free(qm);
  queue_free(qo);
  }

/* -------------------------------------------------------------------------- */
