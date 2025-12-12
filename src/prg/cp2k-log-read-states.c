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
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read one set of electronic states from cp2k output file (MO coefficients)
 
   d        - cp2k data structure
   s        - cp2k state data structure
   n_states - number of states
   f        - open file stream */
void cp2k_log_read_states_set_mo(struct cp2k_dat *d, struct cp2k_state **s, 
  unsigned n_states, FILE *f) {
  unsigned i,j,k,l,i1,i2,n_blocks,n_cols,n_vals = 4;
  char *line;
  struct cp2k_state *state;
  /* initialization */
  n_blocks = (n_states%n_vals ? n_states/n_vals+1 : n_states/n_vals);
  if (!(*s))
    (*s) = cp2k_state_new(n_states);
  state = (*s);
  for (i=0; i<n_states; i++)
    if (!state[i].ao_pj)
      state[i].ao_pj = vec_falloc(d->n_ibfce);
  /* read data blocks */
  for (i=0; i<n_blocks; i++) {
    n_cols = ((i+1)<n_blocks ? n_vals :
      (n_states%n_vals ? n_states%n_vals : n_vals));
    i1 = 0;
    /* header */
    for (j=0; j<2; j++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      line = str_free(line);
      }
    /* energies */
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cp2k log file",1);
    i2 = n_vals*i;
    for (k=0; k<n_cols; k++)
      if (!sscanf(line+23+k*13,"%lf",&(state[i2++].energy)))
        msg_error("invalid format of electronic-state data block",1);
    line = str_free(line);
    /* occupancies */
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cp2k log file",1);
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cp2k log file",1);
    i2 = n_vals*i;
    for (k=0; k<n_cols; k++)
      if (!sscanf(line+23+k*13,"%lf",&(state[i2++].occup)))
        msg_error("invalid format of electronic-state data block",1);
    line = str_free(line);
    /* atoms */
    for (j=0; j<d->n_atoms; j++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error("unexpected end of cp2k log file",1);
      line = str_free(line);
      /* molecular orbitals */
      for (k=0; k<d->kind[d->atom[j].kind].n_sphr_fce; k++) {
        line = str_read_line_new(f);
        if (!line)
          msg_error("unexpected end of cp2k log file",1);
        /* MO coefficients */
        i2 = n_vals*i;
        for (l=0; l<n_cols; l++)
          if (!sscanf(line+23+l*13,"%lf",&(state[i2++].ao_pj[i1])))
            msg_error("invalid format of electronic-state data block",1);
        line = str_free(line);
        i1++;
        }
      }
    }
  }

/* read one set of electronic states from cp2k output file (energies)
 
   d        - cp2k data structure
   s        - cp2k state data structure
   n_states - number of states
   f        - open file stream */
void cp2k_log_read_states_set_en(struct cp2k_dat *d, struct cp2k_state **s, 
  unsigned n_states, FILE *f) {
  unsigned i;
  char *line;
  struct cp2k_state *state;
  /* memory allocation */
  if (!(*s))
    (*s) = cp2k_state_new(n_states);
  state = (*s);
  /* header */
  for (i=0; i<2; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cp2k log file",1);
    line = str_free(line);
    }
  /* energies */
  for (i=0; i<n_states; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of cp2k log file",1);
    if (sscanf(line+11,"%lf%lf",&(state[i].energy),&(state[i].occup))!=2)
      msg_error("invalid format of electronic-state data block",1);
    line = str_free(line);
    }
  }

/* read electronic states from cp2k output file
 
   d - cp2k data struct
   f - open file stream */
void cp2k_log_read_states(struct cp2k_dat *d, FILE *f) {
  unsigned i,j;
  char *line;
  short kpt = 1;
  rewind(f);
  /* sanity check */
  if (d->n_sphr_fce && d->n_ibfce!=d->n_sphr_fce)
    msg_error_f("number of independent AOs (%d) is not equal to"
     " number of spherical functions (%d)",1,d->n_ibfce,d->n_sphr_fce);
  /* k-points */
  if (!d->n_kpoints) {
    d->n_kpoints = 1;
    kpt = 0;
    }
  if (!d->kpoint) {
    d->kpoint = cp2k_kpoint_new(d->n_kpoints);
    if (!kpt) {
      d->kpoint->weight = 1.0;
      for (i=0; i<3; i++)
        d->kpoint->crd[i] = 0.0;
      }
    }
  for (i=0; i<d->n_kpoints; i++) {
    /* memory allocation */
    d->kpoint[i].state = vec_talloc(sizeof(struct cp2k_state*),d->n_spins);
    for (j=0; j<d->n_spins; j++)
      d->kpoint[i].state[j] = NULL;
    /* spin-polarized data */
    if (d->n_spins>1) {
      line = str_ffind_b_new(f," ALPHA MO EIGENVALUES");
      if (line) {
        /* KS energies */
        if (kpt)
          cp2k_log_read_states_set_en(d,
            &(d->kpoint[i].state[0]),d->n_states[0],f);
        /* MO coefficients */
        else
          cp2k_log_read_states_set_mo(d,
            &(d->kpoint[i].state[0]),d->n_states[0],f);
        }
      line = str_free(line);
      line = str_ffind_b_new(f," BETA MO EIGENVALUES");
      if (line) {
        /* KS energies */
        if (kpt)
          cp2k_log_read_states_set_en(d,
            &(d->kpoint[i].state[1]),d->n_states[1],f);
        /* MO coefficients */
        else
          cp2k_log_read_states_set_mo(d,
            &(d->kpoint[i].state[1]),d->n_states[1],f);
        }
      line = str_free(line);
      }
    /* no spin */
    else {
      line = str_ffind_b_new(f," MO EIGENVALUES");
      if (line) {
        /* KS energies */
        if (kpt)
          cp2k_log_read_states_set_en(d,
            &(d->kpoint[i].state[0]),d->n_states[0],f);
        /* MO coefficients */
        else
          cp2k_log_read_states_set_mo(d,
            &(d->kpoint[i].state[0]),d->n_states[0],f);
        }
      line = str_free(line);
      }
    }
  }

/* read electronic states from cp2k binary data file
 
   d    - cp2k data struct
   name - name of the file */
void cp2k_log_read_states_bin(struct cp2k_dat *d, char *name) {
  long unsigned f_size,rec_len;
  int n,n_atoms,n_mos,n_recs;
  unsigned i,j,k;
  double *e,*o,*c;
  FILE *f;
  /* size of the file */
  f_size = file_size(name);
  /* open the file */
  f = file_open(name,"rb");
  /* first-record data */
  file_bin_frt_read_int(&n_atoms,1,f);
  file_bin_frt_read_int(&n_mos,1,f);
  /* sanity check */
  if (d->n_sphr_fce && d->n_ibfce!=d->n_sphr_fce)
    msg_error_f("number of independent AOs (%d) is not equal to"
     " number of spherical functions (%d)",1,d->n_ibfce,d->n_sphr_fce);
  /* data-record length */
  rec_len = 0;
  rec_len += /* number of atoms */
    3*FILE_BIN_FRT_SIZE_INT;
  rec_len += /* number of MOs */
    3*FILE_BIN_FRT_SIZE_INT;
  rec_len += /* MO energies */
    n_mos*FILE_BIN_FRT_SIZE_REAL+2*FILE_BIN_FRT_SIZE_INT;
  rec_len += /* occupation numbers */
    n_mos*FILE_BIN_FRT_SIZE_REAL+2*FILE_BIN_FRT_SIZE_INT;
  rec_len += /* MO coefficients */
    d->n_ibfce*(n_mos*FILE_BIN_FRT_SIZE_REAL+2*FILE_BIN_FRT_SIZE_INT);
  n_recs = f_size/rec_len;
  if (f_size%rec_len)
    msg_error_f("invalid length of the file (file size = %ld,"
      " record length = %ld)",1,f_size,rec_len);
  if (n_recs<d->n_spins || n_recs%d->n_spins)
    msg_error_f("number of data records (%d) incompatible with"
      " number of spin components (%d)",1,n_recs,d->n_spins);
  /* move pointer to data */
  if (fseek(f,f_size-d->n_spins*rec_len,SEEK_SET)<0)
    msg_error("cannot set data reading pointer",1);
  /* sanity check */
  for (i=0; i<d->n_spins; i++)
    if (d->n_states[i]!=n_mos)
      msg_error_f("number of spin-%d states (%d) incompatible with"
        " number of MOs (%d)",1,i+1,d->n_states[i],n_mos);
  /* auxiliary arrays */
  e = vec_falloc(n_mos);
  o = vec_falloc(n_mos);
  c = vec_falloc(n_mos);
  /* k-points */
  if (!d->n_kpoints)
    d->n_kpoints = 1;
  if (d->n_kpoints>1)
    msg_error("reading binary MO coefficients for K-points"
     " is not supported yet",1);
  if (!d->kpoint)
    d->kpoint = cp2k_kpoint_new(d->n_kpoints);
  for (i=0; i<d->n_kpoints; i++) {
    /* memory allocation */
    d->kpoint[i].state = vec_talloc(sizeof(struct cp2k_state*),d->n_spins);
    for (j=0; j<d->n_spins; j++) {
      d->kpoint[i].state[j] = cp2k_state_new(d->n_states[j]);
      for (k=0; k<d->n_states[j]; k++)
        d->kpoint[i].state[j][k].ao_pj = vec_falloc(d->n_ibfce);
      }
    /* alpha spin */
    file_bin_frt_read_int(&n,1,f);
    if (n!=n_atoms)
      msg_error("incompatible number of atoms in alpha-spin data record",1);
    file_bin_frt_read_int(&n,1,f);
    if (n!=n_mos)
      msg_error("incompatible number of MOs in alpha-spin data record",1);
    file_bin_frt_read_real(e,n_mos,f);
    file_bin_frt_read_real(o,n_mos,f);
    for (j=0; j<n_mos; j++) {
      d->kpoint[i].state[0][j].energy = e[j];
      d->kpoint[i].state[0][j].occup = o[j];
      }
    for (j=0; j<d->n_ibfce; j++) {
      file_bin_frt_read_real(c,n_mos,f);
      for (k=0; k<n_mos; k++)
        d->kpoint[i].state[0][k].ao_pj[j] = c[k];
      }
    /* beta spin */
    if (d->n_spins>1) {
      file_bin_frt_read_int(&n,1,f);
      if (n!=n_atoms)
        msg_error("incompatible number of atoms in beta-spin data record",1);
      file_bin_frt_read_int(&n,1,f);
      if (n!=n_mos)
        msg_error("incompatible number of MOs in beta-spin data record",1);
      file_bin_frt_read_real(e,n_mos,f);
      file_bin_frt_read_real(o,n_mos,f);
      for (j=0; j<n_mos; j++) {
        d->kpoint[i].state[1][j].energy = e[j];
        d->kpoint[i].state[1][j].occup = o[j];
        }
      for (j=0; j<d->n_ibfce; j++) {
        file_bin_frt_read_real(c,n_mos,f);
        for (k=0; k<n_mos; k++)
          d->kpoint[i].state[1][k].ao_pj[j] = c[k];
        }
      }
    }
  /* close the file */
  file_close(f);
  /* clean memory */
  vec_ffree(e);
  vec_ffree(o);
  vec_ffree(c);
  }

/* -------------------------------------------------------------------------- */
