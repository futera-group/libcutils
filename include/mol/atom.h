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

#ifndef ZF_LIB_MOL_ATOM_H
#define ZF_LIB_MOL_ATOM_H

/* -------------------------------------------------------------------------- */

#define ATOM_MAX_NUM 111   /* largest valid atomic number */

/* -------------------------------------------------------------------------- */

/* convert atomic number to atom name */
char *atom_name(unsigned);
/* convert atomic number to atom full name */
char *atom_name_full(unsigned);
/* convert PDB atomic symbol to atom name */
char* atom_name_pdb(char*);
/* convert atomic name to PDB position style */
char* atom_name_to_pdb(char*);

/* convert atom name to atomic number */
unsigned atom_num(char*);
/* convert atom name in PDB format to atomic number */
unsigned atom_num_pdb(char*);
/* convert atomic mass to atomic number */
unsigned atom_num_mass(double);

/* convert atomic number to atomic mass */
double atom_mass(unsigned);
/* convert atom name to atomic mass */
double atom_mass_name(char*);
/* convert PDB atom name to atomic mass */
double atom_mass_name_pdb(char*);

/* convert atomic number to atomic radius */
double atom_radius(unsigned);

/* -------------------------------------------------------------------------- */

#endif
