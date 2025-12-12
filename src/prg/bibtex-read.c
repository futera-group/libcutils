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

#include <dirent.h>
#include <stdlib.h>
#include <cmn/dir.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* split string of names to individual authors
 
   list  - the string of names
   item  - name of item (author,editor,...)
   entry - name of entry (label)
   file  - name of bibtex file
   n     - number of names (output) */
char*** bibtex_read_item_author(char *list, char *item , char *entry,
  char *file, unsigned *n) {
  char **av,***s = NULL;
  short error = 0;
  unsigned i,bn,an;
  /* split record to individual authors */
  av = str_split_w(list," and ",&an);
  /* split author to family and given name */
  s = vec_talloc(sizeof(char**),an);
  for (i=0; i<an; i++) {
    s[i] = str_split(av[i],',',&bn);
    if (bn!=2) {
      if (item && entry && file)
        msg_error_f("invalid name \"%s\" in \"%s\" item of \"%s\" entry in"
          " \"%s\" file",1,av[i],item,entry,file);
      else {
        error = 1;
        break;
        }
      }
    str_trim(s[i][0]);
    str_trim(s[i][1]);
    }
  /* clean memory */
  vec_sfree(av,an);
  /* finish */
  if (error) {
    for (i=0; i<an; i++)
      vec_sfree(s[i],2);
    an = 0;
    s = NULL;
    }
  (*n) = an;
  return(s);
  }

/* split string with file data specification into single data categories
 
   bt    - pointer to bibtex database data struct
   rec   - record with file name, path and type specification
   file  - name of bibtex file */
void bibtex_read_item_file(struct bib_entry *bt, char *rec, char *file) {
  unsigned i,sn,nf;
  char **sv,**sf;
  sf = str_split(rec,';',&nf);
  bt->file->n_files = nf;
  bt->file->path = vec_talloc(nf,sizeof(char*));
  bt->file->name = vec_talloc(nf,sizeof(char*));
  bt->file->type = vec_sialloc(nf);
  for (i=0; i<nf; i++) {
    sv = str_split(sf[i],':',&sn);
    if (sn!=2)
      msg_error_f("invalid file specification in \"%s\" entry of \"%s\" file",1,
        bt->name,file);
    bt->file->type[i] = bibtex_type_file_id(sv[1]);
    if (bt->file->type[i]==BIB_FILE_UNKNOWN)
      msg_warn_f("unknown file type \"%s\" in \"%s\" entry of \"%s\" file",
        sv[1],bt->name,file);
    bt->file->path[i] = str_file_dir_new(sv[0]);
    bt->file->name[i] = str_file_name_new(sv[0]);
    sv = vec_sfree(sv,sn);
    }
  sf = vec_sfree(sf,nf);
  }

/* read one item from bibtex entry and save data
 
   bt    - pointer to the bibtex entry data struct
   type  - type of item (author,year,pages,...)
   data  - value of the item */
void bibtex_read_item_assign(struct bib_entry *bt, char *type, char* data) {
  char frmt[]="invalid \"%s\" item value \"%s\" of \"%s\" entry in \"%s\" file";
  unsigned i,n;
  char **v;
  switch (bibtex_type_item_id(type)) {
    case BIB_ITEM_AUTHOR:
      bt->author = bibtex_read_item_author(data,type,
        bt->name,bt->bibfile,&(bt->n_authors));
      break;
    case BIB_ITEM_TITLE:
      bt->title = str_copy_new(data);
      break;
    case BIB_ITEM_JOURNAL:
      bt->journal = str_copy_new(data);
      break;
    case BIB_ITEM_YEAR:
      if (sscanf(data,"%u",&(bt->year))!=1)
        msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      break;
    case BIB_ITEM_PAGES:
      v = str_split(data,'-',&n);
      if (n!=1 && n!=2)
        msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      for (i=0; i<2; i++)
        if (i<n && sscanf(v[i],"%u",&(bt->page[i]))!=1)
           msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      vec_sfree(v,n);
      break;
    case BIB_ITEM_KEYWORD:
      bt->keyword = str_split(data,',',&(bt->n_keywords));
      for (i=0; i<bt->n_keywords; i++) {
        str_trim(bt->keyword[i]);
        }
      break;
    case BIB_ITEM_OWNER:
      bt->owner = str_copy_new(data);
      break;
    case BIB_ITEM_TIMESTAMP:
      if (sscanf(data,"%u.%u.%u",&(bt->timestamp[2]),&(bt->timestamp[1]),
        &(bt->timestamp[0]))!=3)
        msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      break;
    case BIB_ITEM_VOLUME:
      if (sscanf(data,"%u",&(bt->volume))!=1)
        msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      break;
    case BIB_ITEM_FILE:
      bibtex_read_item_file(bt,data,bt->bibfile);
      break;
    case BIB_ITEM_NUMBER:
      if (sscanf(data,"%u",&(bt->number))!=1)
        msg_error_f(frmt,1,type,data,bt->name,bt->bibfile);
      break;
    case BIB_ITEM_ORGANIZ:
      bt->organization = str_copy_new(data);
      break;
    case BIB_ITEM_SCHOOL:
      bt->school = str_copy_new(data);
      break;
    case BIB_ITEM_ADDRESS:
      bt->address = str_copy_new(data);
      break;
    case BIB_ITEM_CHAPTER:
      bt->chapter = str_copy_new(data);
      break;
    case BIB_ITEM_PUBLISH:
      bt->publisher = str_copy_new(data);
      break;
    case BIB_ITEM_EDITOR:
      bt->editor = bibtex_read_item_author(data,type,
        bt->name,bt->bibfile,&(bt->n_editors));
      break;
    case BIB_ITEM_BOOKTITLE:
      bt->book_title = str_copy_new(data);
      break;
    case BIB_ITEM_UNKNOWN:
      msg_warn_f("unknown item \"%s\" of \"%s\" entry in \"%s\" file",
        type,bt->name,bt->bibfile);
      break;
    }
  }

/* read one item from bibtex entry and save data
 
   bt    - pointer to the bibtex entry data struct
   data  - line with item data
   file  - pointer to open database file */
void bibtex_read_item(struct bib_entry *bt, char *data, FILE *file) {
  char *type = NULL,*drec = NULL,*line = NULL,*s = NULL,*t = NULL;
  unsigned in,id,nd,nl;
  short ok = 0;
  /* type of item */
  in = 0;
  id = 0;
  nd = str_length(data);
  type = str_new(nd);
  while (id<nd && data[id]!='=' && data[id]!=' ' && 
         data[id]!='\t' && data[id]!='\n')
    type[in++] = data[id++];
  type[in]='\0';
  /* data record */
  in = 0;
  id++;
  drec = str_new(nd);
  while (id<nd && data[id]!='{')
    id++;
  id++;
  while (id<nd && data[id]!='}')
    drec[in++] = data[id++];
  drec[in] = '\0';
  if (id && id<nd && data[id]=='}')
    ok = 1;
  /* multiple-line data record */
  while (!ok) {
    line = str_read_line_new(file);
    if (!line)
      msg_error_f("unexpected end of file while reading"
        " \"%s\" item of \"%s\" entry in \"%s\" file",1,
        type,bt->name,bt->bibfile);
    str_trim(line);
    id = 0;
    nl = str_length(line);
    s = str_new(nl);
    while (id<nl && line[id]!='}') {
      s[id] = line[id];
      id++;
      }
    s[id] = '\0';
    t = str_new(str_length(drec)+str_length(s)+1);
    sprintf(t,"%s %s",drec,s);
    str_free(drec);
    drec = t;
    if (id && id<nd && line[id]=='}')
      ok = 1;
    str_free(line);
    str_free(s);
    }
  /* save data */
  bibtex_read_item_assign(bt,type,drec);
  /* clean memory */
  str_free(type);
  str_free(drec);
  }

/* read one entry from BibTeX database file 
 
   db    - bibtex database
   head  - type of the entry (first line)
   fname - file name
   file  - pointer to open database file */
void bibtex_read_entry(struct list *db, char *head, char *fname, FILE *file) {
  char *line = NULL,*t = NULL;
  unsigned ih,in,nh;
  struct bib_entry *bt;
  bt = bibtex_entry_new();
  bt->bibfile = str_copy_new(fname);
  /* type of entry */
  in = 0;
  ih = 1;
  nh = str_length(head);
  t = str_new(nh);
  while (ih<nh && head[ih]!='{' && head[ih]!=' ' && 
         head[ih]!='\t' && head[ih]!='\n')
    t[in++] = head[ih++];
  t[in] = '\0';
  bt->type = bibtex_type_entry_id(t);
  if (bt->type!=BIB_TYPE_COMMENT && bt->type!=BIB_TYPE_UNKNOWN) {
    /* name of entry */
    in = 0;
    ih++;
    while (ih<nh && head[ih]!=',' && head[ih]!=' ' &&
           head[ih]!='\t' && head[ih]!='\n')
      t[in++] = head[ih++];
    t[in] = '\0';
    bt->name = str_copy_new(t);
    /* entry items */
    for (line=str_read_line_new(file); line;
         line=str_free(line),line=str_read_line_new(file)) {
      str_trim(line);
      if (line[0]=='}') {
        line = str_free(line);
        break;
        }
      bibtex_read_item(bt,line,file);
      }
    /* add entry to database */
    list_add_end_p(db,bt);
    }
  else {
    if (bt->type==BIB_TYPE_UNKNOWN)
      msg_warn_f("unknown type \"%s\" of entry \"%s\" in \"%s\" file",
        t,bt->name,fname);
    bibtex_entry_free(bt);
    }
  /* clean memory */
  str_free(t);
  }

/* read all entries from one BibTeX database file
 
   db   - bibtex database
   name - name of the file */
void bibtex_read_file(struct list *db, char *name) {
  char *line;
  FILE *file;
  file = file_open(name,"r");
  for (line=str_read_line_new(file); line;
       line=str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (line[0]=='@')
      bibtex_read_entry(db,line,name,file);
    }
  file_close(file);
  }

/* read all BibTeX database files in specified directory
  
   path - path to directory with bibtex files
   name - name of the bibtex file 
   nf   - number of bibtex files (output) */
struct list* bibtex_read_dir(char *path, char *name, unsigned *nf) {
  char *full = NULL,*type = NULL;
  DIR *dir;
  struct dirent *dir_info;
  struct stat file_stat;
  struct list *db;
  (*nf) = 0;
  db = list_alloc();
  /* one specified file */
  if (name && name[0]) {
    bibtex_read_file(db,name);
    (*nf) = 1;
    }
  /* search directory */
  else {
    if (path && path[0])
      dir = dir_open(path,1);
    else
      dir = dir_open_current();
    dir_info = readdir(dir);
    while (dir_info) {
      if (path && path[0]) {
        full = str_new(str_length(path)+str_length(dir_info->d_name)+1);
        sprintf(full,"%s/%s",path,dir_info->d_name);
        }
      else
        full = str_copy_new(dir_info->d_name);
      lstat(full,&file_stat);
      if (S_ISREG(file_stat.st_mode) || S_ISLNK(file_stat.st_mode)) {
        type = str_file_ext_new(dir_info->d_name);
        if (str_compare(type,"bib")) {
          bibtex_read_file(db,full);
          (*nf)++;
          }
        type = str_free(type);
        }
      dir_info = readdir(dir);
      full = str_free(full);
      }
    dir_close(dir);
    }
  return(db);
  }

/* -------------------------------------------------------------------------- */
