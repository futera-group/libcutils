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

#include <cmn/file.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/units.h>
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read electronic couplings from cp2k binary data file
 
   d    - cp2k data struct
   name - name of the file
   what - what to read (0 = all, n = one block, -1 = no couplings */
void cp2k_log_read_ecoupl_bin(struct cp2k_dat *d, char *name, short what) {
  int n_spins,ni,*vi;
  unsigned i,j,k,l,n,id;
  double *vf;
  FILE *f;
  /* open the file */
  f = file_open(name,"rb");
  /* number of spin components */
  file_bin_frt_read_int(&n_spins,1,f);
  if (d->n_spins!=n_spins)
    msg_error_f("inconsistent number of spin components in log/hab files"
      " (%d/%d)",1,d->n_spins,n_spins);
  /* memory allocation */
  if (!d->ecpl)
    d->ecpl = cp2k_ecpl_new();
  if (!d->ecpl->n_states) {
    d->ecpl->n_states = vec_ualloc(d->n_spins);
    vec_uset(d->ecpl->n_states,0,d->n_spins);
    }
  if (!d->ecpl->energy) {
    d->ecpl->energy = mat_falloc(d->n_spins,0);
    mat_fset(d->ecpl->energy,0.0,d->n_spins,0);
    }
  /* number of POD blocks */
  file_bin_frt_read_int(&ni,1,f);
  if (d->ecpl->n_blocks) {
    if (d->ecpl->n_blocks != ni)
      msg_error_f("inconsisten number of POD blocks in LOG and HAB files"
        " (%d vs %d)",1,d->ecpl->n_blocks,ni);
    }
  else {
    d->ecpl->n_blocks = ni;
    d->ecpl->block = cp2k_ecpl_block_new(d->ecpl->n_blocks);
    }
  /* block dimensions */
  vi = vec_ialloc(d->ecpl->n_blocks);
  file_bin_frt_read_int(vi,d->ecpl->n_blocks,f);
  for (i=0; i<d->ecpl->n_blocks; i++) {
    if (!what || what<0 || what==(i+1)) {
      /* number of states */
      if (!d->ecpl->block[i].n_states) {
        d->ecpl->block[i].n_states = vec_ualloc(d->n_spins);
        for (j=0; j<d->n_spins; j++)
          d->ecpl->block[i].n_states[j] = vi[i];
        }
      else {
        for (j=0; j<d->n_spins; j++)
          if (d->ecpl->block[i].n_states[j] != vi[i])
            msg_error_f("inconsistent number of POD states in block %d"
              " in LOG and HAB files (%d vs %d)",1,i+1,
              d->ecpl->block[i].n_states[j],vi[i]);
        }
      /* state energies */
      if (!d->ecpl->block[i].energy) {
        d->ecpl->block[i].energy = vec_talloc(sizeof(double*),d->n_spins);
        for (j=0; j<d->n_spins; j++) {
          d->ecpl->block[i].energy[j] = NULL;
          }
        }
      }
    /* coupling elements */
    if (!what || what==(i+1)) {
      if (!d->ecpl->block[i].coupling) {
        d->ecpl->block[i].coupling = vec_talloc(sizeof(double***),d->n_spins);
        for (j=0; j<d->n_spins; j++)
          d->ecpl->block[i].coupling[j] = NULL;
        }
      }
    }
  vi = vec_ifree(vi);
  /* spin components */
  for (i=0; i<d->n_spins; i++) {
    /* state energies */
    for (j=0; j<d->ecpl->n_blocks; j++) {
      if (!what || what<0 || what==(j+1)) {
        if (!d->ecpl->block[j].energy[i])
          d->ecpl->block[j].energy[i] = 
            vec_falloc(d->ecpl->block[j].n_states[i]);
        file_bin_frt_read_real(d->ecpl->block[j].energy[i],
          d->ecpl->block[j].n_states[i],f);
        }
      else {
        vf = vec_falloc(d->ecpl->block[j].n_states[i]);
        file_bin_frt_read_real(vf,d->ecpl->block[j].n_states[i],f);
        vf = vec_ffree(vf);
        }
      }
    /* coupling elements */
    for (j=0; j<d->ecpl->n_blocks; j++) {
      if (!what || what==(j+1)) {
        if (!d->ecpl->block[j].coupling[i])
          d->ecpl->block[j].coupling[i] = 
            vec_talloc(sizeof(double**),d->ecpl->n_blocks-j-1);
        for (k=0; k<d->ecpl->n_blocks-j-1; k++)
          d->ecpl->block[j].coupling[i][k] = NULL;
        }
      vf = vec_falloc(d->ecpl->block[j].n_states[i]);
      for (k=j+1; k<d->ecpl->n_blocks; k++) {
        id = k-j-1;
        if ((!what || what==(j+1)) && !d->ecpl->block[j].coupling[i][id]) {
          d->ecpl->block[j].coupling[i][id] = 
            mat_falloc(d->ecpl->block[j].n_states[i],
              d->ecpl->block[k].n_states[i]);
          }
        for (l=0; l<d->ecpl->block[k].n_states[i]; l++) {
          file_bin_frt_read_real(vf,d->ecpl->block[j].n_states[i],f);
          if (!what || what==(j+1)) {
            for (n=0; n<d->ecpl->block[j].n_states[i]; n++) {
              d->ecpl->block[j].coupling[i][id][n][l] = vf[n];
              }
            }
          }
        }
      vec_ffree(vf);
      }
    }
  /* close the file */
  file_close(f);
  }

/* read localized POD wavefunctions from cp2k binary data file
 
   d    - cp2k data struct
   name - name of the file
   what - what to read (0 = all, n = one block, -1 = no coefficients */
void cp2k_log_read_ecoupl_wfn(struct cp2k_dat *d, char *name, short what) {
  int n_spins,ni,*vi;
  unsigned i,j,k;
  double *vf;
  FILE *f;
  /* open the file */
  f = file_open(name,"rb");
  /* number of spin components */
  file_bin_frt_read_int(&n_spins,1,f);
  if (d->n_spins!=n_spins)
    msg_error_f("inconsistent number of spin components in log/pod-wfn files"
      " (%d/%d)",1,d->n_spins,n_spins);
  /* memory allocation */
  if (!d->ecpl)
    d->ecpl = cp2k_ecpl_new();
  /* number of POD blocks */
  file_bin_frt_read_int(&ni,1,f);
  if (!d->ecpl->n_blocks) {
    d->ecpl->n_blocks = ni;
    d->ecpl->block = cp2k_ecpl_block_new(d->ecpl->n_blocks);
    }
  else if (d->ecpl->n_blocks!=ni)
    msg_error_f("inconsistent number of POD blocks (%d/%d)",1,
      d->ecpl->n_blocks,ni);
  /* block dimensions */
  vi = vec_ialloc(d->ecpl->n_blocks);
  file_bin_frt_read_int(vi,d->ecpl->n_blocks,f);
  for (i=0; i<d->ecpl->n_blocks; i++) {
    if (!d->ecpl->block[i].n_states) {
      d->ecpl->block[i].n_states = vec_ualloc(d->n_spins);
      for (j=0; j<d->n_spins; j++)
        d->ecpl->block[i].n_states[j] = vi[i];
      }
    else {
      for (j=0; j<d->n_spins; j++)
        if (d->ecpl->block[i].n_states[j]!=vi[i])
          msg_error_f("inconsistent number of POD block-%d/spin-%d states"
            " (%d/%d)",1,i+1,j+1,d->ecpl->block[i].n_states[j],vi[i]);
      }
    }
  vi = vec_ifree(vi);
  /* wavefunctions */
  if (what>=0) {
    for (i=0; i<d->ecpl->n_blocks; i++) 
      if ((!what || what==(i+1)) && !d->ecpl->block[i].ao_pj) {
        d->ecpl->block[i].ao_pj = vec_talloc(sizeof(double**),d->n_spins);
        for (j=0; j<d->n_spins; j++)
          d->ecpl->block[i].ao_pj[j] = NULL;
        }
    for (i=0; i<d->n_spins; i++) {
      for (j=0; j<d->ecpl->n_blocks; j++) {
        /* save data */
        if (!what || what==(j+1)) {
          if (!d->ecpl->block[j].ao_pj[i])
            d->ecpl->block[j].ao_pj[i] = 
              mat_falloc(d->ecpl->block[j].n_states[i],d->n_ibfce);
          for (k=0; k<d->ecpl->block[j].n_states[i]; k++) 
            file_bin_frt_read_real(d->ecpl->block[j].ao_pj[i][k],d->n_ibfce,f);
          }
        /* skip data */
        else {
          vf = vec_falloc(d->n_ibfce);
          for (k=0; k<d->ecpl->block[j].n_states[i]; k++)
            file_bin_frt_read_real(vf,d->n_ibfce,f);
          vf = vec_ffree(vf);
          }
        }
      }
    }
  /* close the file */
  file_close(f);
  }

/* read POD transformation matrix from cp2k binary data file
 
   d    - cp2k data struct
   name - name of the file
   what - type of the matrix (1 = forward, 2 = reverse) */
void cp2k_log_read_ecoupl_tmx(struct cp2k_dat *d, char *name, short what) {
  unsigned i,j;
  double *vf;
  FILE *f;
  /* open the file */
  f = file_open(name,"rb");
  /* memory allocation */
  d->tmat = mat_falloc(d->n_ibfce,d->n_ibfce);
  /* read matrix elements */
  vf = vec_falloc(d->n_ibfce);
  for (i=0; i<d->n_ibfce; i++) {
    file_bin_frt_read_real(vf,d->n_ibfce,f);
    for (j=0; j<d->n_ibfce; j++)
      d->tmat[j][i] = vf[j];
    }
  vf = vec_ffree(vf);
  /* close the file */
  file_close(f);
  }

/* read electronic-coupling block states from cp2k output file
 
   d        - cp2k data struct
   n_states - number of states
   energy   - state energies
   f        - open file stream */
void cp2k_log_read_ecoupl_states(struct cp2k_dat *d, unsigned **n_states,
  double ***energy, FILE *f) {
  char *line = NULL; 
  unsigned i,n;
  double oc1,oc2,en1,en2;
  struct queue *q1,*q2;
  /* header */
  line = str_read_line_new(f);
  line = str_free(line);
  /* read state energies */
  q1 = queue_alloc();
  q2 = queue_alloc();
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    if (str_length(line)<10) {
      line = str_free(line);
      break;
      }
    /* spin-polarized calculation */
    if (d->n_spins>1) {
      if (sscanf(line,"%u%lf%lf%lf%lf",&n,&oc1,&en1,&oc2,&en2)==5) {
        queue_fadd(q1,en1);
        queue_fadd(q2,en2);
        }
      else if (sscanf(line,"%u%lf%lf",&n,&oc1,&en1)!=3)
        msg_error("invalid format of total-system state energies",1);
      }
    /* spin-restricted calculation */
    else {
      if (sscanf(line,"%u%lf%lf",&n,&oc1,&en1)!=3)
        msg_error("invalid format of total-system state energies",1);
      queue_fadd(q1,en1);
      }
    }
    if (!(*n_states))
      (*n_states) = vec_ualloc(d->n_spins);
    if (!(*energy)) {
      (*energy) = vec_talloc(sizeof(double*),d->n_spins);
      for (i=0; i<d->n_spins; i++)
        (*energy)[i] = NULL;
      }
    (*n_states)[0] = q1->num;
    if (!((*energy)[0]))
      (*energy)[0] = vec_falloc((*n_states)[0]);
    for (i=0; i<(*n_states)[0]; i++)
      queue_fget(q1,&((*energy)[0][i]));
    if (d->n_spins>1) {
      (*n_states)[1] = q2->num;
      if (!((*energy)[1]))
        (*energy)[1] = vec_falloc((*n_states)[1]);
      for (i=0; i<(*n_states)[1]; i++)
        queue_fget(q2,&((*energy)[1][i]));
      }
  /* clean memory */
  queue_free(q1);
  queue_free(q2);
  }

/* read electronic couplings from cp2k output file
 
   d - cp2k data struct
   w - what to read
   f - open file stream */
void cp2k_log_read_ecoupl(struct cp2k_dat *d, int w, FILE *f) {
  char *line,flag[80];
  unsigned i,j,k,n,m,*v;
  unsigned ib1,ib2,is1,is2,ix,id;
  double cp1,cp2;
  struct queue *q1;
  rewind(f);
  /* find the electronic-state data block */
  line = str_ffind_new(f,"Electronic coupling - Projection-operator method");
  if (line) {
    line = str_free(line);
    /* memory allocation */
    d->ecpl = cp2k_ecpl_new();
    /* number of blocks */
    line = str_ffind_new(f,"Number of fragments");
    if (!line) 
      msg_error("number of coupling blocks is not specified",1);
    if (str_length(line)<39 || sscanf(line+38,"%u",&(d->ecpl->n_blocks))!=1)
      msg_error("cannot read number of coupling blocks",1);
    line = str_free(line);
    /* memory allocation */
    d->ecpl->block = cp2k_ecpl_block_new(d->ecpl->n_blocks);
    /* block data */
    for (i=0; i<d->ecpl->n_blocks; i++) {
      sprintf(flag,"  Block %d:",i+1);
      line = str_ffind_b_new(f,flag);
      if (!line)
        msg_error_f("cannot found coupling block #%d specification",1,i+1);
      line = str_free(line);
      /* number of atoms */
      line = str_read_line_new(f);
      if (!line || !str_sub_bfind(line,"  Number of block atoms") ||
          str_length(line)<39 || 
          sscanf(line+38,"%u",&(d->ecpl->block[i].n_atoms))!=1)
        msg_error_f("cannot read number of coupling block #%d atoms",1,i+1);
      line = str_free(line);
      /* number of electrons */
      line = str_read_line_new(f);
      if (!line || !str_sub_bfind(line,"  Number of block electrons") ||
          str_length(line)<39 || 
          sscanf(line+38,"%u",&(d->ecpl->block[i].n_electrons))!=1)
        msg_error_f("cannot read number of coupling block #%d electrons",1,i+1);
      /* atom IDs */
      for (j=0; j<2; j++) {
        line = str_free(line);
        line = str_read_line_new(f);
        }
      if (!line || !str_sub_bfind(line,"  Block atom IDs"))
        msg_error_f("cannot find atom IDs in coupling block #%d",1,i+1);
      /* in-line IDs */
      if (d->ecpl->block[i].n_atoms<10) {
        if (str_length(line)<39)
          msg_error("invalid format of coupling block IDs",1);
        d->ecpl->block[i].atom_id = str_parse_uarray(line+38,&n);
        if (d->ecpl->block[i].n_atoms!=n)
          msg_error_f("inconsistent number of atom IDs"
            " in coupling block #%d (%d/%d)",1,i+1,
            d->ecpl->block[i].n_atoms,n);
        line = str_free(line);
        }
      /* multiple-line IDs */
      else {
        q1 = queue_alloc();
        n = d->ecpl->block[i].n_atoms/10;
        if (d->ecpl->block[i].n_atoms%10)
          n++;
        line = str_free(line);
        for (j=0; j<n; j++) {
          line = str_read_line_new(f);
          v = str_parse_uarray(line,&m);
          for (k=0; k<m; k++)
            queue_uadd(q1,v[k]);
          vec_ufree(v);
          line = str_free(line);
          }
        if (d->ecpl->block[i].n_atoms!=q1->num)
          msg_error_f("inconsistent number of atom IDs"
            " in coupling block #%d (%d/%d)",1,i+1,
             d->ecpl->block[i].n_atoms,q1->num);
        d->ecpl->block[i].atom_id = vec_ualloc(d->ecpl->block[i].n_atoms);
        for (j=0; j<d->ecpl->block[i].n_atoms; j++) {
          queue_uget(q1,&(d->ecpl->block[i].atom_id[j]));
          d->ecpl->block[i].atom_id[j]--;
          }
        queue_free(q1);
        }
      }
    if (w & READ_ECOUPL) {
      /* system state energies */
      line = str_ffind_b_new(f,"  State energies (the whole system):");
      if (!line) 
        msg_error("cannot read state energies",1);
      line = str_free(line);
      cp2k_log_read_ecoupl_states(d,&(d->ecpl->n_states),&(d->ecpl->energy),f);
      /* block state energies */
      for (i=0; i<d->ecpl->n_blocks; i++) {
        sprintf(flag,"  State energies (block %d states):",i+1);
        line = str_ffind_b_new(f,flag);
        if (!line) 
          msg_error_f("cannot read coupling block #%d state energies",1,i+1);
        line = str_free(line);
        cp2k_log_read_ecoupl_states(d,&(d->ecpl->block[i].n_states),
          &(d->ecpl->block[i].energy),f);
        }
      /* couplings */
      line = str_ffind_b_new(f,"  Coupling elements [meV]:");
      if (line) {
        line = str_free(line);
        line = str_read_line_new(f);
        /* memory allocation */
        for (ib1=0; ib1<d->ecpl->n_blocks; ib1++) {
          d->ecpl->block[ib1].coupling = vec_talloc(sizeof(double***),d->n_spins);
          for (ix=0; ix<d->n_spins; ix++) {
            d->ecpl->block[ib1].coupling[ix] = 
              vec_talloc(sizeof(double**),d->ecpl->n_blocks-ib1-1);
            for (ib2=ib1+1; ib2<d->ecpl->n_blocks; ib2++) {
              id = ib2-ib1-1;
              d->ecpl->block[ib1].coupling[ix][id] = 
                vec_talloc(sizeof(double*),d->ecpl->block[ib1].n_states[ix]);
              for (is1=0; is1<d->ecpl->block[ib1].n_states[ix]; is1++) {
                d->ecpl->block[ib1].coupling[ix][id][is1] = 
                  vec_falloc(d->ecpl->block[id].n_states[ix]);
                vec_fset(d->ecpl->block[ib1].coupling[ix][id][is1],
                  0.0,d->ecpl->block[id].n_states[ix]);
                }
              }
            }
          }
        /* read coupling values */
        for (line=str_free(line),line=str_read_line_new(f); line;
             line=str_free(line),line=str_read_line_new(f)) {
          if (str_length(line)<20) {
            line = str_free(line);
            break;
            }
          /* elements IDs */
          if (sscanf(line,"%d[%d] - %d[%d]",&ib1,&is1,&ib2,&is2)!=4)
            msg_error("invalid format of electronic coupling elements",1);
          if (ib1<1 || ib1>d->ecpl->n_blocks)
            msg_error_f("invalid coupling block-1 ID (%d)\n",ib1);
          ib1--;
          if (ib2<1 || ib2>d->ecpl->n_blocks)
            msg_error_f("invalid coupling block-2 ID (%d)\n",ib2);
          ib2--;
          if (is1<1 || is2<1)
            msg_error("invalid coupling state ID\n",1);
          is1--;
          is2--;
          id = ib2-ib1-1;
          /* elements value */
          if (sscanf(line+22,"%lf%lf",&cp1,&cp2)==2) {
            if (is1<d->ecpl->block[ib1].n_states[0] && 
                is2<d->ecpl->block[ib2].n_states[0])
              d->ecpl->block[ib1].coupling[0][id][is1][is2] = 
                cp1/CONV_H_EV/1000.0;
            if (is1<d->ecpl->block[ib1].n_states[1] && 
                is2<d->ecpl->block[ib2].n_states[1])
              d->ecpl->block[ib1].coupling[1][id][is1][is2] = 
                cp2/CONV_H_EV/1000.0;
            }
          else if (sscanf(line+22,"%lf",&cp1)==1) {
            if (is1<d->ecpl->block[ib1].n_states[0] && 
                is2<d->ecpl->block[ib2].n_states[0])
              d->ecpl->block[ib1].coupling[0][id][is1][is2] = 
                cp1/CONV_H_EV/1000.0;
            }
          else
            msg_error_f("cannot read coupling element %d[%d]-%d[%d]",
              1,ib1,is1,ib2,is2);
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
