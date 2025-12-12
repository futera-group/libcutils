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
#include <cmn/string.h>
#include <cmn/print.h>
#include <cmn/queue.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* print one entry number item from bibtex database
 
   num   - the number
   label - printing label */
void bibtex_print_number(unsigned num, char *label) {
  printf("%-*s",BIB_PRNT_INDT_WDTH,label);
  printf("%u\n",num);
  }

/* print one entry text item from bibtex database
 
   text  - item contents
   label - printing label */
void bibtex_print_text(char *text, char *label) {
  unsigned i,i0 = 0,i1,n,m,nl = 0;
  if (!text)
    return;
  printf("%-*s",BIB_PRNT_INDT_WDTH,label);
  n = str_length(text);
  m = i1 = (BIB_PRNT_LINE_WDTH-BIB_PRNT_INDT_WDTH);
  while (i0<n) {
    if (nl)
      printf("\n%*c",BIB_PRNT_INDT_WDTH,' ');
    if ((n-i0+1)<=m) {
      printf("%s",text+i0);
      break;
      }
    while (i0<n && text[i0]==' ')
      i0++;
    while (i1>1 && text[i1]!=' ')
      i1--;
    if (i1>(i0+(m/2)))
      for (i=i0; i<i1; i++)
        printf("%c",text[i]);
    else {
      printf("%-*.*s",m,m,text+i0);
      i1 = i0+m;
      }
    i0 = i1+1;
    i1 += m;
    nl++;
    }
  printf("\n");
  }

/* print pages from one bibtex database entry
 
   page  - array with page numbers */
void bibtex_print_pages(unsigned *page) {
  if (page && page[0]) {
    printf("%-*s",BIB_PRNT_INDT_WDTH,"Pages:");
    printf("%d",page[0]);
    if (page[1])
      printf("-%d",page[1]);
    printf("\n");
    }
  }

/* print one entry authors from bibtex database
 
   name    - name array
   n_names - number of names
   label   - printing label */
void bibtex_print_author(char ***name, unsigned n_names, char *label) {
  char ss[10],*st = NULL,*sn = NULL,*sl = NULL;
  unsigned i,k,n;
  if (name) {
    /* create string */
    for (i=0; i<n_names; i++) {
      /* initials of given names */
      k = 0;
      n = str_length(name[i][1]);
      while (k<n) {
        while (k<n && name[i][1][k]==' ')
          k++;
        if (k<n && name[i][1][k]>='A' && name[i][1][k]<='Z') {
          sprintf(ss,"%c.",name[i][1][k]);
          st = str_merge_new(sn,ss);
          str_free(sn);
          sn = st;
          }
        while (k<n && name[i][1][k]!=' ')
          k++;
        }
      /* surname */
      st = str_merge_new(sn,name[i][0]);
      str_free(sn);
      sn = st;
      /* comma */
      if (i>0) {
        sprintf(ss,"%s",", ");
        st = str_merge_new(ss,sn);
        str_free(sn);
        sn = st;
        }
      /* name list */
      st = str_merge_new(sl,sn);
      sn = str_free(sn);
      str_free(sl);
      sl = st;
      }
    /* print to screen */
    bibtex_print_text(sl,label);
    str_free(sl);
    }
  }

/* print one entry citation from bibtex database
 
   jrnl - name of journal
   vol  - journal volume
   page - pages
   year - year */
void bibtex_print_ref(char *jrnl, unsigned vol, unsigned *page, unsigned year) {
  char st[80],*s = NULL,*t = NULL;
  unsigned i,n;
  struct queue *q;
  q = queue_alloc();
  /* journal */
  if (jrnl && jrnl[0]) {
    n = str_length(jrnl);
    for (i=0; i<n; i++)
      if (jrnl[i]!=' ' && !(jrnl[i]=='\\' && i<(n-1) && jrnl[i+1]=='&')) {
        queue_cadd(q,jrnl[i]);
        }
    if (q->num) {
      s = str_new(q->num);
      n = 0;
      while (q->num)
        queue_cget(q,&(s[n++]));
      s[n] = '\0';
      }
    }
  /* volume */
  if (vol) {
    if (s)
      sprintf(st," %d",vol);
    else
      sprintf(st,"%d",vol);
    t = str_merge_new(s,st);
    str_free(s);
    s = t;
    }
  /* page */
  if (page && page[0]) {
    if (s) {
      if (page[1])
        sprintf(st," %d-%d",page[0],page[1]);
      else
        sprintf(st," %d",page[0]);
      }
    else {
      if (page[1])
        sprintf(st,"%d-%d",page[0],page[1]);
      else
        sprintf(st,"%d",page[0]);
      }
    t = str_merge_new(s,st);
    str_free(s);
    s = t;
    }
  /* year */
  if (year) {
    if (s)
      sprintf(st," (%d)",year);
    else
      sprintf(st,"(%d)",year);
    t = str_merge_new(s,st);
    str_free(s);
    s = t;
    }
  /* print out reference */
  bibtex_print_text(s,"Journal:");
  }

/* print one entry keywords from bibtex database

   key    - keyword array
   n_keys - number of keywords
   label  - printing label */
void bibtex_print_key(char **key, unsigned n_keys, char *label) {
  char *s = NULL,*t = NULL;
  unsigned i;
  if (n_keys) {
    for (i=0; i<n_keys; i++) {
      if (i>0) {
        t = str_merge_new(s,", ");
        str_free(s);
        s = t;
        }
      t = str_merge_new(s,key[i]);
      str_free(s);
      s = t;
      }
    bibtex_print_text(s,label);
    }
  }

/* print info recorded in one bibtex database entry 
 
   e - pointer to the entry
   m - print out list of missing fields
   f - print the entry in frame */
void bibtex_print_entry(struct bib_entry *e, short m, short f) {
  /* frame */
  if (f)
    print_Hline('*',80);
  /* bibtex entry */
  bibtex_print_text(bibtex_type_entry_name(e->type),"Type:");
  bibtex_print_author(e->author,e->n_authors,"Authors:");
  bibtex_print_author(e->editor,e->n_editors,"Editors:");
  if (e->title && e->title[0])
    bibtex_print_text(e->title,"Title:");
  if (e->book_title && e->book_title[0])
    bibtex_print_text(e->book_title,"Book title:");
  if (e->chapter && e->chapter[0])
    bibtex_print_text(e->chapter,"Chapter:");
  if (e->journal && e->journal[0])
    bibtex_print_text(e->journal,"Journal:");
  if (e->publisher && e->publisher[0])
    bibtex_print_text(e->publisher,"Publisher:");
  if (e->organization && e->organization[0])
    bibtex_print_text(e->organization,"Organiz.:");
  if (e->school && e->school[0])
    bibtex_print_text(e->school,"School:");
  if (e->address && e->address[0])
    bibtex_print_text(e->address,"Address:");
  if (e->year)
    bibtex_print_number(e->year,"Year:");
  if (e->volume)
    bibtex_print_number(e->volume,"Volume:");
  if (e->number)
    bibtex_print_number(e->number,"Number:");
  bibtex_print_pages(e->page);
  bibtex_print_key(e->keyword,e->n_keywords,"Keywords:");
  /* frame */
  if (f)
    print_Hline('-',80);
  /* missign fields */
  if (m)
    bibtex_entry_check(e);
  /* frame */
  if (f)
    print_Hline('*',80);
  }

/* -------------------------------------------------------------------------- */
