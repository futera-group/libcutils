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
#include "cmn/file.h"
#include "cmn/matrix.h"
#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"
#include "cmn/types.h"
#include "cmn/units.h"
#include "cmn/vector.h"
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* read keyword and return its internal ID

   line - keyword line from the data file
   key  - extracted keyword (output) 
   opt  - optional keyword argument value (output) */
short mld_read_key(char *line, char *key, char *opt) {
  char **s;
  short id;
  unsigned n;
  s = str_split(line,' ',&n);
  if (n<1 || str_length(s[0])<3 ||
      s[0][0] != '[' || s[0][str_length(s[0])-1] != ']')
    msg_error_f("invalid data format: %s\n",1,line);
  /* keyword */
  str_copy(key,s[0]+1);
  key[str_length(key)-1] = '\0';
  id = mld_data_type_id(key);
  /* optional arguments */
  if (id == MLD_DATA_ATOM) {
    if (n < 2) {
      msg_warn("coordinate units not specified, Angstroms assumed");
      sprintf(opt,"AU");
      }
    else
      str_copy(opt,s[1]);
    }
  else if (id == MLD_DATA_GEOC) {
    if (n < 2) {
      msg_warn("geometry format not specified, XYZ assumed");
      sprintf(opt,"XYZ");
      }
    else
      str_copy(opt,s[1]);
    }
  vec_sfree(s,n);
  return(id);
  }

/* read atomic coordinates from the molden data file

   m    - molden data struct
   unit - units ("Angs" or "AU")
   f    - open file stream */
char* mld_read_atoms(struct mld_dat *m, char *unit, FILE *f) {
  char *line,**s;
  unsigned i,j,n;
  short conv = 0;
  struct mld_atom *a;
  struct queue *q;
  q = queue_alloc();
  /* read atomic data to list */
  for (line = str_read_line_new(f); line;
       line = str_free(line), line = str_read_line_new(f)) {
     s = str_split(line,' ',&n);
     if (n < 1)
       continue;
     if (s[0][0] == '[')
       break;
     if (n < 6)
       msg_error("invalid format of atomic coordinates",1);
     a = mld_atom_new(1);
     a->sym = str_copy_new(s[0]);
     if (!type_read(s[2],&(a->num),TYPE_UINT))
       msg_error("invalid format of atomic coordinates",1);
     if (!type_read(s[3],&(a->crd[0]),TYPE_DOUBLE))
       msg_error("invalid format of atomic coordinates",1);
     if (!type_read(s[4],&(a->crd[1]),TYPE_DOUBLE))
       msg_error("invalid format of atomic coordinates",1);
     if (!type_read(s[5],&(a->crd[2]),TYPE_DOUBLE))
       msg_error("invalid format of atomic coordinates",1);
     queue_add(q,a);
     vec_sfree(s,n);
     }
  /* unit conversion */
  if (str_compare(unit,"AU"))
    conv = 0;
  else if (str_compare(unit,"Angs"))
    conv = 1;
  else
    msg_error_f("unknown coordinate units \"%s\"",1,unit);
  /* save data to the molden structure */
  m->n_atoms = q->num;
  m->atom = mld_atom_new(m->n_atoms);
  for (i=0; i<m->n_atoms; i++) {
    a = queue_get(q);
    m->atom[i].sym = a->sym;
    a->sym = NULL;  
    m->atom[i].num = a->num;
    for (j=0; j<3; j++) {
      m->atom[i].crd[j] = a->crd[j];
      if (conv)
        m->atom[i].crd[j] /= CONV_B_ANG;
      }
    }
  queue_free(q);
  return(line);
  }

/* read vibrational frequencies from the molden data file

   m - molden data struct
   f - open file stream */
char* mld_read_freq(struct mld_dat *m, FILE *f) {
  char *line,**s;
  unsigned i,n;
  double v;
  struct queue *q;
  q = queue_alloc();
  /* read atomic data to list */
  for (line = str_read_line_new(f); line;
       line = str_free(line), line = str_read_line_new(f)) {
    s = str_split(line,' ',&n);
    if (n < 1)
      continue;
    if (s[0][0] == '[')
      break;
    if (!type_read(s[0],&v,TYPE_DOUBLE))
      msg_error("invalid format of vibrational frequencies",1);
    queue_fadd(q,v);
    vec_sfree(s,n);
    }
  /* save data to the molden structure */
  m->n_vib_modes = q->num;
  m->freq = mld_freq_new(m->n_vib_modes);
  for (i=0; i<m->n_vib_modes; i++)
    queue_fget(q,&(m->freq[i].freq));
  queue_free(q);
  return(line);
  }

/* read vibrational normal-mode atomic dispacements from the molden data file

   m - molden data struct
   f - open file stream */
char* mld_read_freq_modes(struct mld_dat *m, FILE *f) {
  char *line = NULL,**s;
  unsigned i,j,n,id;
  for (i=0; i<m->n_vib_modes; i++) {
    /* mode ID */
    line = str_free(line);
    line = str_read_line_new(f);
    s = str_split(line,' ',&n);
    if (n < 2 || !str_compare(s[0],"vibration") ||
        !type_read(s[1],&id,TYPE_UINT) || id != i+1)
      msg_error("invalid format of vibrational normal modes",1);
    /* atomic displacements */
    m->freq[i].displ = mat_falloc(m->n_atoms,3);
    for (j=0; j<m->n_atoms; j++) {
      line = str_free(line);
      line = str_read_line_new(f);
      if (sscanf(line,"%lf%lf%lf",
          &(m->freq[i].displ[j][0]),
          &(m->freq[i].displ[j][1]),
          &(m->freq[i].displ[j][2])) != 3)
        msg_error("invalid format of vibrational normal modes",1);
      }
    }
  return(line);
  }

/* read intensities of vibrational frequencies from the molden data file

   m - molden data struct
   f - open file stream */
char* mld_read_freq_ints(struct mld_dat *m, FILE *f) {
  char *line,**s;
  unsigned i,n;
  double v;
  struct queue *q;
  q = queue_alloc();
  /* read atomic data to list */
  for (line = str_read_line_new(f); line;
       line = str_free(line), line = str_read_line_new(f)) {
    s = str_split(line,' ',&n);
    if (n < 1)
      continue;
    if (s[0][0] == '[')
      break;
    if (!type_read(s[0],&v,TYPE_DOUBLE))
      msg_error("invalid format of vibrational frequency intensities",1);
    queue_fadd(q,v);
    vec_sfree(s,n);
    }
  /* save data to the molden structure */
  if (q->num != m->n_vib_modes)
    msg_error("invalid number of vibrational frequency intentisites",1);
  for (i=0; i<m->n_vib_modes; i++)
    queue_fget(q,&(m->freq[i].ir));
  queue_free(q);
  return(line);
  }

/* read data from molden data file
 
   m - pointer to molden data struct
   f - name of the log file */
void mld_read(struct mld_dat *m, char *f) {
  char *line,key[80],opt[80];
  FILE *file;
  /* open the file */
  file = file_open(f,"r");
  /* find beginning of the data record */
  line = str_ffind_new(file,"[Molden Format]");
  if (!line)
    msg_error_f("no data record found in \"%s\" file",1,f);
  /* read the file */
  line = str_free(line);
  line = str_read_line_new(file);
  while (line) {
    str_trim(line);
    /* keyword */
    if (str_length(line) && line[0]=='[') {
      switch (mld_read_key(line,key,opt)) {
        case 0: 
          msg_error_f("unknown data type \"%s\"",1,key); break;
        case MLD_DATA_ATOM:
          line = mld_read_atoms(m,opt,file); continue;
        case MLD_DATA_FREQ:
          line = mld_read_freq(m,file); continue;
        case MLD_DATA_FRNC:
          line = mld_read_freq_modes(m,file); continue;
        case MLD_DATA_FRI:
          line = mld_read_freq_ints(m,file); continue;
        default:
          msg_warn_f("skipping data \"%s\"",key); break;
        }
      }
    line = str_free(line);
    line = str_read_line_new(file);
    }
  /* close the file */
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
