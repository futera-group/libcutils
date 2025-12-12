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

#ifndef ZF_LIB_PRG_BIBTEX_H
#define ZF_LIB_PRG_BIBTEX_H

#include <stdio.h>
#include <cmn/list.h>

/* symbolic constants */
#define BIB_MAX_ITEM          21   /* maximal number of items */

/* bibtex entry types */
#define BIB_TYPE_COMMENT       0   /* comment */
#define BIB_TYPE_ARTICLE       1   /* article */
#define BIB_TYPE_BOOK          2   /* book */
#define BIB_TYPE_BOOKLET       3   /* booklet */
#define BIB_TYPE_CONF          4   /* conference */
#define BIB_TYPE_INBOOK        5   /* part of book */
#define BIB_TYPE_INCOLL        6   /* part of collection */
#define BIB_TYPE_INPROC        7   /* part of proceedings */
#define BIB_TYPE_MANUAL        8   /* manual */
#define BIB_TYPE_MTHESIS       9   /* master thesis */
#define BIB_TYPE_MISC         10   /* others */
#define BIB_TYPE_PHDTHESIS    11   /* phd thesis */
#define BIB_TYPE_PROC         12   /* proceedings */
#define BIB_TYPE_TECHREP      13   /* technical report */
#define BIB_TYPE_UNPUB        14   /* unpublished */
#define BIB_TYPE_UNKNOWN      15   /* unknown type */

/* bibtex item types */
#define BIB_ITEM_LABEL         0   /* label */
#define BIB_ITEM_TYPE          1   /* label */
#define BIB_ITEM_AUTHOR        2   /* author */
#define BIB_ITEM_TITLE         3   /* title */
#define BIB_ITEM_JOURNAL       4   /* journal */
#define BIB_ITEM_YEAR          5   /* year */
#define BIB_ITEM_PAGES         6   /* pages */
#define BIB_ITEM_OWNER         7   /* owner */
#define BIB_ITEM_TIMESTAMP     8   /* timestamp */
#define BIB_ITEM_VOLUME        9   /* volume */
#define BIB_ITEM_FILE         10   /* file */
#define BIB_ITEM_NUMBER       11   /* number */
#define BIB_ITEM_ORGANIZ      12   /* organization */
#define BIB_ITEM_SCHOOL       13   /* school */
#define BIB_ITEM_ADDRESS      14   /* address */
#define BIB_ITEM_CHAPTER      15   /* chapter */
#define BIB_ITEM_PUBLISH      16   /* publisher */
#define BIB_ITEM_EDITOR       17   /* editor */
#define BIB_ITEM_BOOKTITLE    18   /* book title */
#define BIB_ITEM_KEYWORD      19   /* keyword */
#define BIB_ITEM_UNKNOWN      20   /* unknown item */

/* bibtex file types */
#define BIB_FILE_PDF           1   /* PDF */
#define BIB_FILE_UNKNOWN       2   /* unknown file type */

/* bibtex journal name format */
#define BIB_JRNL_NAME_FILE     1   /* file name style */

/* printing */
#define BIB_PRNT_INDT_WDTH    12   /* indent width */
#define BIB_PRNT_LINE_WDTH    80   /* line width */

/* bibtex database entry - file data struct */
struct bib_entry_file {
  unsigned n_files;                /* number of files */
  char **path;                     /* path to the file */
  char **name;                     /* name of the file */
  short *type;                     /* type of the file */
  };

/* bibtex database entry data struct */
struct bib_entry {
  short type;                      /* type of entry (article,report,...) */
  char ***author;                  /* name of authors */
  unsigned n_authors;              /* number of authors */
  char ***editor;                  /* name of editors */
  unsigned n_editors;              /* number of editors */
  char **keyword;                  /* keywords */
  unsigned n_keywords;             /* number of keywords */
  char *name;                      /* name of entry (label) */
  char *title;                     /* title */
  char *book_title;                /* book title */
  char *chapter;                   /* book chapter */
  char *journal;                   /* journal */
  char *publisher;                 /* publisher */
  char *organization;              /* organization */
  char *school;                    /* school */
  char *owner;                     /* owner */
  char *address;                   /* address */
  unsigned year;                   /* year */
  unsigned volume;                 /* volume */
  unsigned number;                 /* number */
  unsigned page[2];                /* pages xx-yy */
  unsigned timestamp[3];           /* time stamp dd/mm/yy */
  char *bibfile;                   /* name of bibtex database file */
  struct bib_entry_file *file;     /* file (pdf,tex,txt,...) */
  };

/* free memory allocated for list of authors */
char*** bibtex_entry_author_free(char***, unsigned n);

/* allocate memory for new bibtex entry file data structure */
struct bib_entry_file* bibtex_entry_file_new(void);
/* free memory allocated for bibtex entry file data structure */
struct bib_entry_file* bibtex_entry_file_free(struct bib_entry_file*);

/* allocate memory for new bibtex database entry data structure */
struct bib_entry* bibtex_entry_new(void);
/* free memory allocated for bibtex database entry data structure */
struct bib_entry* bibtex_entry_free(struct bib_entry*);

/* print info recorded in one bibtex database entry */
void bibtex_print_entry(struct bib_entry*, short, short);
/* print one entry keywords from bibtex database */
void bibtex_print_key(char**, unsigned, char*);
/* print one entry citation from bibtex database */
void bibtex_print_ref(char*, unsigned, unsigned*, unsigned);
/* print one entry authors from bibtex database */
void bibtex_print_author(char***, unsigned, char*);
/* print one entry text item from bibtex database */
void bibtex_print_text(char*, char*);

/* return internal ID of specified bibtex entry type */
short bibtex_type_entry_id(char*);
/* return name of speciefied bibtex entry type */
char* bibtex_type_entry_name(short);
/* return internal ID of specified bibtex item type */
short bibtex_type_item_id(char*);
/* return internal ID of specified bibtex file type */
short bibtex_type_file_id(char*type);
/* return name of specified bibtex file type */
char* bibtex_type_file_name(short);
/* return name of specified bibtex entry item */
char* bibtex_type_item_name(short);

/* return string with bibtex entry name generated from the saved data */
char* bibtex_entry_name_new(struct bib_entry*);
/* check mandatory options and print out missing ones */
void bibtex_entry_check(struct bib_entry*);

/* read one item from bibtex entry and save data */
void bibtex_read_item_assign(struct bib_entry*, char*, char*);
/* split string of names to individual authors */
char*** bibtex_read_item_author(char*, char*, char*, char*, unsigned*);
/* split string with file data specification into single data categories */
void bibtex_read_item_file(struct bib_entry*, char*, char*);
/* read one item from bibtex entry and save data */
void bibtex_read_item(struct bib_entry*, char*, FILE*);
/* read one entry from BibTeX database file */
void bibtex_read_entry(struct list*, char*, char*, FILE*);
/* read all entries from one BibTeX database file */
void bibtex_read_file(struct list*, char*);
/* read all BibTeX database files in specified directory */
struct list* bibtex_read_dir(char*, char*, unsigned*);

/* write list of bibtex entries to one file */
void bibtex_write(struct list*, char*, short);
/* write one bibtex entry to file */
void bibtex_write_entry(struct bib_entry*, FILE*);
/* write one bibtex entry to file */
void bibtex_write_entry_f(struct bib_entry*, char*, short);

#endif
