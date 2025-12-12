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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* read atom specification from open apc file
 
   p - pointer to amber prep struct
   f - pointer to open amber prep file */
void mol_apc_fread_atoms(struct apc_mol *p, FILE *f) {
  char *line,atom_name[80],atom_type[80],tree_mark[10];
  unsigned id = 0;
  struct queue *q;
  struct apc_atom *at;
  /* read data to queue */
  q = queue_alloc();
  for (line=str_read_line_new(f); line;
       str_free(line),line=str_read_line_new(f)) {
    at = mol_apc_atom_new(1);
    if (sscanf(line,"%u%s%s%s%lf%lf%lf%lf",&(at->id),atom_name,atom_type,
      tree_mark,&(at->coord[0]),&(at->coord[1]),&(at->coord[2]),
      &(at->charge))!=8) {
      mol_apc_atom_free(at,1);
      break;
      }
    at->name = str_copy_new(atom_name);
    at->type = str_copy_new(atom_type);
    at->tree = mol_apc_tree_id(tree_mark);
    queue_add(q,at);
    }
  /* convert queue to array */
  p->n_atoms = q->num;
  p->atom = mol_apc_atom_new(p->n_atoms);
  while (q->num) {
    at = queue_get(q);
    mol_apc_atom_copy(at,p->atom+id);
    mol_apc_atom_free(at,1);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read loop specification from open apc file
 
   p - pointer to amber prep struct
   f - pointer to open amber prep file */
void mol_apc_fread_loops(struct apc_mol *p, FILE *f) {
  char *line,s1[80],s2[80];
  unsigned *b,id = 0;
  struct queue *q;
  /* find data block */
  line = str_ffind_new(f,"LOOP");
  if (!line)
    return;
  /* read data to queue */
  q = queue_alloc();
  for (line=str_read_line_new(f); line;
       str_free(line),line=str_read_line_new(f)) {
    if (sscanf(line,"%s%s",s1,s2)!=2)
      break;
    b = vec_ualloc(2);
    b[0] = mol_apc_atom_id(p,s1);
    b[1] = mol_apc_atom_id(p,s2);
    queue_add(q,b);
    }
  /* convert queue to array */
  p->n_loops = q->num;
  p->loop = mat_ualloc(p->n_loops,2);
  while (q->num) {
    b = queue_get(q);
    p->loop[id][0] = b[0];
    p->loop[id][1] = b[1];
    vec_ufree(b);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read improper specification from open apc file
 
   p - pointer to amber prep struct
   f - pointer to open amber prep file */
void mol_apc_fread_imprs(struct apc_mol *p, FILE *f) {
  char *line,s[4][80];
  unsigned *b,i,id = 0;
  struct queue *q;
  /* find data block */
  line = str_ffind_new(f,"IMPROPER");
  if (!line)
    return;
  /* read data to queue */
  q = queue_alloc();
  for (line=str_read_line_new(f); line;
       str_free(line),line=str_read_line_new(f)) {
    if (sscanf(line,"%s%s%s%s",s[0],s[1],s[2],s[3])!=4)
      break;
    b = vec_ualloc(4);
    for (i=0; i<4; i++)
      b[i] = mol_apc_atom_id(p,s[i]);
    queue_add(q,b);
    }
  /* convert queue to array */
  p->n_imprs = q->num;
  p->impr = mat_ualloc(p->n_imprs,4);
  while (q->num) {
    b = queue_get(q);
    for (i=0; i<4; i++)
      p->impr[id][i] = b[i];
    vec_ufree(b);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read molecular data from open apc file
 
   p - pointer to amber prep struct
   f - pointer to open amber prep file */
void mol_apc_fread(struct apc_mol *p, FILE *f) {
  char *line,s1[80],s2[80];
  unsigned i;
  /* comment */
  line = str_read_line_new(f);
  if (!line)
    msg_error("amber prep file is empty",1);
  str_free(line);
  /* database name */
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of prep file while reading database name",1);
  str_free(line);
  /* title */
  p->title = str_read_line_new(f);
  if (!p->title)
    msg_error("unexpected end of prep file while reading title",1);
  str_trim(p->title);
  /* file name */
  p->file = str_read_line_new(f);
  if (!p->file)
    msg_error("unexpected end of prep file while reading residuum file name",1);
  str_trim(p->file);
  /* residuum name */
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of prep file while reading residuum name",1);
  if (sscanf(line,"%s%s",s1,s2)!=2)
    msg_error("invalid format of prep file - residuum name expected",1);
  p->resname = str_copy_new(s1);
  line = str_free(line);
  /* coordinate format */
  if (!str_compare(s2,"XYZ"))
    msg_error("cartesian coordinate type prep file format expected",1);
  /* skip header */
  for (i=0; i<5; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of prep file while reading header",1);
    line = str_free(line);
    }
  /* atom specification */
  mol_apc_fread_atoms(p,f);
  /* loop specification */
  mol_apc_fread_loops(p,f);
  /* improper dihedrals */
  mol_apc_fread_imprs(p,f);
  }

/* read molecular data from apc file
 
   p - pointer to amber prep struct
   f - name of amber prep file */
void mol_apc_read(struct apc_mol *p, char *f) {
  FILE *file = file_open(f,"r");
  mol_apc_fread(p,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */

/* write molecular data to file in apc format
 
   p     - amber prep data structure
   first - first-residuum database indicator
   last  - last-residuum database indicator
   f     - open file stream for output */
void mol_apc_fwrite(struct apc_mol *p, short first, short last, FILE *f) {
  short sgn[3][3] = {{1,1,-1},{1,-1,1},{-1,1,1}};
  char tree_mark[10],s[4][80];
  unsigned i,j;
  /* header */
  if (first)
    fprintf(f,"    0    0    2\n\n");
  fprintf(f,"%s\n",p->title);
  fprintf(f,"%s\n",p->file);
  fprintf(f,"%-3s   XYZ  0\n",p->resname);
  fprintf(f,"CHANGE     OMIT DU   BEG\n");
  fprintf(f,"%8.4f\n",0.0);
  for (i=0; i<3; i++)
    fprintf(f,"%4d  DUMM  DU    M%15.3f%12.3f%12.3f%15.3f\n",
      i+1,999.0*sgn[i][0],999.0*sgn[i][1],999.0*sgn[i][2],0.0);
  /* atom specification */
  for (i=0; i<p->n_atoms; i++) {
    mol_apc_tree_mark(p->atom[i].tree,tree_mark);
    fprintf(f,"%4d  %-4s  %-6s%-4s%15.6f%12.6f%12.6f%12.6f\n",
      p->atom[i].id,p->atom[i].name,p->atom[i].type,tree_mark,
      p->atom[i].coord[0],p->atom[i].coord[1],p->atom[i].coord[2],
      p->atom[i].charge);
    }
  fprintf(f,"\n");
  /* loop specification */
  if (p->n_loops) {
    fprintf(f,"\nLOOP\n");
    for (i=0; i<p->n_loops; i++) {
      for (j=0; j<2; j++)
        str_trim_copy(p->atom[p->loop[i][j]].name,s[j]);
      fprintf(f,"%5s%5s\n",s[0],s[1]);
      }
    fprintf(f,"\n");
    }
  /* improper specification */
  if (p->n_imprs) {
    fprintf(f,"\nIMPROPER\n");
    for (i=0; i<p->n_imprs; i++) {
      for (j=0; j<4; j++)
        str_trim_copy(p->atom[p->impr[i][j]].name,s[j]);
      fprintf(f,"%5s%5s%5s%5s\n",s[0],s[1],s[2],s[3]);
      }
    fprintf(f,"\n");
    }
  fprintf(f,"DONE\n");
  if (last)
    fprintf(f,"STOP\n");
  }

/* write molecular data to file in apc format
 
   p     - amber prep data structure
   first - first-residuum database indicator
   last  - last-residuum database indicator
   f     - name of amber prep file */
void mol_apc_write(struct apc_mol *p, short first, short last, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  mol_apc_fwrite(p,first,last,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
