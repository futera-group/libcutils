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

#ifndef ZF_LIB_MOL_AACID_H
#define ZF_LIB_MOL_AACID_H

/* -------------------------------------------------------------------------- */

/* amino acids */
#define AA_ALA   1   /* alanine */
#define AA_ARG   2   /* arginine */
#define AA_ASN   3   /* asparagine */
#define AA_ASP   4   /* aspartic acid */
#define AA_CYS   5   /* cysteine */
#define AA_GLU   6   /* glutamic acid */
#define AA_GLN   7   /* glutamine */
#define AA_GLY   8   /* glycine */
#define AA_HIS   9   /* histidine */
#define AA_ILE  10   /* isoleucine */
#define AA_LEU  11   /* leucine */
#define AA_LYS  12   /* lysine */
#define AA_MET  13   /* methionine */
#define AA_PHE  14   /* phenylalanine */
#define AA_PRO  15   /* proline */
#define AA_SER  16   /* serine */
#define AA_THR  17   /* threonine */
#define AA_TRP  18   /* tryptophan */
#define AA_TYR  19   /* tyrosine */
#define AA_VAL  20   /* valine */

/* structure forms */
#define AA_FORM_NEUTRAL     0    /* neutral form (NH2,COOH) */
#define AA_FORM_NTER_N      1    /* neutral N-terminus (NH3+) */
#define AA_FORM_CTER_N      2    /* neutral C-terminus (COO-) */
#define AA_FORM_CHAIN_N     3    /* neutral peptide chain */
#define AA_FORM_ZWITTERION  4    /* zwitterion (NH3+,COO-) */
#define AA_FORM_NTER_P      5    /* protonated N-terminus (NH3+) */
#define AA_FORM_CTER_P      6    /* deprotonated C-terminus (COO-) */
#define AA_FORM_CHAIN_P     7    /* de/protonated peptide chain */

/* -------------------------------------------------------------------------- */

#include <cmn/list.h>

/* Convert one-letter amino-acid name to internal ID */
short aacid_id_1(char);
/* Convert three-letter amino-acid name to internal ID */
short aacid_id_3(char*);
/* Convert full amino-acid name to internal ID */
short aacid_id_full(char*);

/* Return 1-letter symbol of specified amino acid */
char* aacid_name_1(short);
/* Return 3-letter name of specified amino acid */
char* aacid_name_3(short);
/* Return full name of specified amino acid */
char* aacid_name_full(short);

/* Create Z-matrix list definition of alanine */
struct list* aacid_geom_def_ala(short);
/* Create Z-matrix list definition of arginine */
struct list* aacid_geom_def_arg(short);
/* Create Z-matrix list definition of asparagine */
struct list* aacid_geom_def_asn(short);
/* Create Z-matrix list definition of aspartic acid */
struct list* aacid_geom_def_asp(short);
/* Create Z-matrix list definition of cysteine */
struct list* aacid_geom_def_cys(short);
/* Create Z-matrix list definition of glutamic acid */
struct list* aacid_geom_def_glu(short);
/* Create Z-matrix list definition of glutamine */
struct list* aacid_geom_def_gln(short);
/* Create Z-matrix list definition of glycine */
struct list* aacid_geom_def_gly(short);
/* Create Z-matrix list definition of histidine */
struct list* aacid_geom_def_his(short);
/* Create Z-matrix list definition of isoleucine */
struct list* aacid_geom_def_ile(short);
/* Create Z-matrix list definition of leucine */
struct list* aacid_geom_def_leu(short);
/* Create Z-matrix list definition of lysine */
struct list* aacid_geom_def_lys(short);
/* Create Z-matrix list definition of methionine */
struct list* aacid_geom_def_met(short);
/* Create Z-matrix list definition phenylalanine */
struct list* aacid_geom_def_phe(short);
/* Create Z-matrix list definition of proline */
struct list* aacid_geom_def_pro(short);
/* Create Z-matrix list definition of serine */
struct list* aacid_geom_def_ser(short);
/* Create Z-matrix list definition of threonine */
struct list* aacid_geom_def_thr(short);
/* Create Z-matrix list definition of tryptophan */
struct list* aacid_geom_def_trp(short);
/* Create Z-matrix list definition of tyrosine */
struct list* aacid_geom_def_tyr(short);
/* Create Z-matrix list definition of valine */
struct list* aacid_geom_def_val(short);

/* Create structure of a single amino acid (list of Z-matrix records) */
struct list* aacid_geom_def(short, short);
/* Create structure of a single amino acid (Z-matrix format) */
struct zmt_mol* aacid_geom_zmt(short, short);

/* Structure modifications (C-terminus,N-terminus,chain) */
void aacid_geom_mod(struct list*, short, short);

/* create chain of amino-acids (peptide structure) */
struct pdb_mol* aacid_chain(char*, short, short);

/* -------------------------------------------------------------------------- */

#endif
