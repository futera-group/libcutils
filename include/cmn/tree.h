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

#ifndef ZF_LIB_CMN_TREE_H
#define ZF_LIB_CMN_TREE_H

/* -------------------------------------------------------------------------- */

/* structs and types */
struct tdata {
  struct tdata **t_child;  /* pointer array to ascendants */
  struct tdata *t_parent;  /* pointer to descandant */
  unsigned t_num;          /* numbef or ascendants */
  unsigned t_id;           /* internal ID */
  void *t_data;            /* pointer to saved data */
  };

struct tree {
  struct tdata *root;      /* pointer to root */
  unsigned num;            /* total number of nodes */
  };

/* -------------------------------------------------------------------------- */

/* allocate new tree structure */
struct tree *tree_alloc(void);
/* allocate memory for tree node data and copy it */
void *tree_alloc_data(short, void*);
/* destroy allocated tree structure */
void tree_free(struct tree*);
/* clean tree structure */
void tree_clean(struct tree*);

/* set root of tree structure */
struct tdata *tree_root_set(struct tree*, short, void*);
/* change root of the tree to specifed node */
void tree_root_change(struct tree*, struct tdata*);
/* reorganize tree to maximalize path from root to leaf */
void tree_root_maxpath(struct tree*);

/* add child node to parent one in tree graph */
struct tdata *tree_child_add(struct tree*, struct tdata*, short, void*);
/* sort childern in the tree node */
void tree_child_sort(struct tdata*, short(f)(void*,void*));
/* increase child vector in tree graph node */
void tree_child_vec_inc(struct tdata*, struct tdata*);
/* delete one child pointer from tree graph node  */
void tree_child_vec_del(struct tdata*, struct tdata*);

/* return length of path from root to given node */
unsigned tree_path_len_p(struct tree*, struct tdata*);
/* return length of longest path in the tree graph */
unsigned tree_path_len_l(struct tree*, struct tdata**);
/* return length of shortest path in the tree graph */
unsigned tree_path_len_s(struct tree*, struct tdata**);

/* sort descendants in tree according to pathway lenghts */
void tree_sort_len(struct tree*);

/* find element in tree according to specified comparing function */
void *tree_find_data_t(struct tree*, void*, short(f)(void*,void*), unsigned*);

/* print out string data saved in the tree structure */
void tree_sprint(struct tree*);

/* -------------------------------------------------------------------------- */

#endif
