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
#include <cmn/list.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* write one bibtex entry to open file
 
   e - the bibtex entry
   f - open stream of the file
   w - rewrite existing file */
void bibtex_write_entry(struct bib_entry *e, FILE *f) {
  unsigned i;
  fprintf(f,"\n");
  /* type and label */
  fprintf(f,"@%s{%s,\n",bibtex_type_entry_name(e->type),e->name);
  /* authors */
  if (e->author && e->n_authors) {
    fprintf(f,"  author = {");
    for (i=0; i<e->n_authors; i++) {
      if (i>0)
        fprintf(f," and ");
      fprintf(f,"%s",e->author[i][0]);
      if (e->author[i][1])
        fprintf(f,", %s",e->author[i][1]);
      }
    fprintf(f,"},\n");
    }
  /* editors */
  if (e->editor && e->n_editors) {
    fprintf(f,"  editor = {");
    for (i=0; i<e->n_editors; i++) {
      if (i>0)
        fprintf(f," and ");
      fprintf(f,"%s",e->editor[i][0]);
      if (e->editor[i][1])
        fprintf(f,", %s",e->editor[i][1]);
      }
    fprintf(f,"},\n");
    }
  /* title */
  if (e->title && e->title[0])
    fprintf(f,"  title = {%s},\n",e->title);
  /* book title */
  if (e->book_title && e->book_title[0])
    fprintf(f,"  booktitle = {%s},\n",e->book_title);
  /* book chapter */
  if (e->chapter && e->chapter[0])
    fprintf(f,"  chapter = {%s},\n",e->chapter);
  /* journal */
  if (e->journal && e->journal[0])
    fprintf(f,"  journal = {%s},\n",e->journal);
  /* publisher */
  if (e->publisher && e->publisher[0])
    fprintf(f,"  publisher = {%s},\n",e->publisher);
  /* organization */
  if (e->organization && e->organization[0])
    fprintf(f,"  organization = {%s},\n",e->organization);
  /* organization */
  if (e->school && e->school[0])
    fprintf(f,"  school = {%s},\n",e->school);
  /* address */
  if (e->address && e->address[0])
    fprintf(f,"  address = {%s},\n",e->address);
  /* year */
  if (e->year)
    fprintf(f,"  year = {%u},\n",e->year);
  /* volume */
  if (e->volume)
    fprintf(f,"  volume = {%u},\n",e->volume);
  /* number */
  if (e->number)
    fprintf(f,"  number = {%u},\n",e->number);
  /* pages */
  if (e->page[0]) {
    fprintf(f,"  pages = {%u",e->page[0]);
    if (e->page[1]) 
      fprintf(f,"-%u",e->page[1]);
    fprintf(f,"},\n");
    }
  /* keyword */
  if (e->keyword && e->n_keywords) {
    fprintf(f,"  keywords = {");
    for (i=0; i<e->n_keywords; i++) {
      if (i>0)
        fprintf(f,", ");
      fprintf(f,"%s",e->keyword[i]);
      }
    fprintf(f,"},\n");
    }
  /* f */
  if (e->file->n_files) {
    fprintf(f,"  file = {");
    for (i=0; i<e->file->n_files; i++) {
      if (i>0)
        fprintf(f,";");
      if (e->file->path && e->file->path[0])
        fprintf(f,":%s/%s:%s",e->file->path[i],e->file->name[i],
          bibtex_type_file_name(e->file->type[i]));
      else
        fprintf(f,":%s:%s",e->file->name[i],
          bibtex_type_file_name(e->file->type[i]));
      }
    fprintf(f,"},\n");
    }
  /* owner */
  if (e->owner && e->owner[0])
    fprintf(f,"  owner = {%s},\n",e->owner);
  /* timestamp */
  if (e->timestamp[0] && e->timestamp[1] && e->timestamp[2])
    fprintf(f,"  timestamp = {%04u.%02u.%02u}\n",
      e->timestamp[2],e->timestamp[1],e->timestamp[0]);
  fprintf(f,"}\n");
  }

/* write one bibtex entry to file
 
   e - the bibtex entry
   f - name of the file
   w - rewrite existing file */
void bibtex_write_entry_f(struct bib_entry *e, char *f, short w) {
  FILE *file;
  file = file_open(f,(w ? "w" : "a"));
  bibtex_write_entry(e,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */

/* write list of bibtex entries to one file
 
   r - the list of bibtex entries
   f - name of the file
   w - rewrite existing file */
void bibtex_write(struct list *r, char *f, short w) {
  struct ldata *p;
  FILE *file;
  /* open file */
  file = file_open(f,(w ? "w" : "a"));
  /* records */
  if (r && r->first)
    for (p=r->first; p; p=p->l_next)
      bibtex_write_entry(p->l_data,file);
  /* close file */
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
