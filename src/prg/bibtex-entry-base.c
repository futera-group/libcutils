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
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for new bibtex entry file data structure */
struct bib_entry_file* bibtex_entry_file_new(void) {
  struct bib_entry_file *b = NULL;
  /* memory allocation */
  b = (struct bib_entry_file*)malloc(sizeof(struct bib_entry_file));
  if (!b)
    msg_error("cannot callote memory for bibtex file data structure",1);
  /* initialization */
  b->n_files = 0;
  b->path = NULL;
  b->name = NULL;
  b->type = NULL;
  return(b);
  }

/* allocate memory for new bibtex database entry data structure */
struct bib_entry* bibtex_entry_new(void) {
  struct bib_entry *b = NULL;
  /* memory allocation */
  b = (struct bib_entry*)malloc(sizeof(struct bib_entry));
  if (!b)
    msg_error("cannot allocate memory for new bibtex database data struct",1);
  /* initialization */
  b->type = 0;
  b->author = NULL;
  b->n_authors = 0;
  b->editor = NULL;
  b->n_editors = 0;
  b->keyword = NULL;
  b->n_keywords = 0;
  b->name = NULL;
  b->title = NULL;
  b->book_title = NULL;
  b->chapter = NULL;
  b->journal = NULL;
  b->publisher = NULL;
  b->organization = NULL;
  b->school = NULL;
  b->owner = NULL;
  b->address = NULL;
  b->year = 0;
  b->volume = 0;
  b->number = 0;
  b->page[0] = 0;
  b->page[1] = 0;
  b->timestamp[0] = 0;
  b->timestamp[1] = 0;
  b->timestamp[2] = 0;
  b->bibfile = NULL;
  b->file = bibtex_entry_file_new();
  return(b);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for list of authors

   a - array of author names
   n - number of authors */
char*** bibtex_entry_author_free(char ***a, unsigned n) {
  unsigned i;
  if (a) {
    for (i=0; i<n; i++)
      vec_sfree(a[i],2);
    free(a);
    a = NULL;
    }
  return(a);
  }

/* free memory allocated for bibtex entry file data structure 

   b - the bibtex entry file */
struct bib_entry_file* bibtex_entry_file_free(struct bib_entry_file *b) {
  if (b) {
    b->path = vec_sfree(b->path,b->n_files);
    b->name = vec_sfree(b->name,b->n_files);
    b->type = vec_sifree(b->type);
    free(b);
    b = NULL;
    }
  return(b);
  }

/* free memory allocated for bibtex database entry data structure

   b - the bibtex entry */
struct bib_entry* bibtex_entry_free(struct bib_entry *b) {
  if (b) {
    b->author = bibtex_entry_author_free(b->author,b->n_authors);
    b->editor = bibtex_entry_author_free(b->editor,b->n_editors);
    if (b->keyword)
      b->keyword = vec_sfree(b->keyword,b->n_keywords);
    b->name = str_free(b->name);
    b->title = str_free(b->title);
    b->book_title = str_free(b->book_title);
    b->chapter = str_free(b->chapter);
    b->journal = str_free(b->journal);
    b->publisher = str_free(b->publisher);
    b->organization = str_free(b->organization);
    b->school = str_free(b->school);
    b->owner = str_free(b->owner);
    b->address = str_free(b->address);
    b->bibfile = str_free(b->bibfile);
    b->file = bibtex_entry_file_free(b->file);
    free(b);
    b = NULL;
    }
  return(b);
  }

/* -------------------------------------------------------------------------- */
