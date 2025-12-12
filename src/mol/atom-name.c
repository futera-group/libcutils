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

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <cmn/string.h>
#include "mol/atom.h"

/* -------------------------------------------------------------------------- */

/* global variables */
char atom_names[ATOM_MAX_NUM][3] = {
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
  "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
  "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg"};

char atom_full_names[ATOM_MAX_NUM][25] = {
  "Hydrogen",     "Helium",      "Lithium",     "Beryllium",     "Boron",
  "Carbon",       "Nitrogen",    "Oxygen",      "Fluorine",      "Neon",
  "Sodium",       "Magnesium",   "Aluminium",   "Silicon",       "Phosphorus",
  "Sulfur",       "Chlorine",    "Argon",       "Potassium",     "Calcium",
  "Scandium",     "Titanium",    "Vanadium",    "Chromium",      "Manganese",
  "Iron",         "Cobalt",      "Nickel",      "Copper",        "Zinc",
  "Gallium",      "Germanium",   "Arsenic",     "Selenium",      "Bromine",
  "Krypton",      "Rubidium",    "Strontium",   "Yttrium",       "Zirconium",
  "Niobium",      "Molybdenium", "Technetium",  "Ruthenium",     "Rhodium",
  "Palladium",    "Silver",      "Cadmium",     "Indium",        "Tin", 
  "Antimony",     "Tellurium",   "Iodine",      "Xenon",         "Cesium", 
  "Barium",       "Lanthanum",   "Cerium",      "Praseodymium",  "Neodymium",
  "Promethium",   "Samarium",    "Europium",    "Gadolinium",    "Terbium", 
  "Dysprosium",   "Holmium",     "Erbium",      "Thulium",       "Ytterbium", 
  "Lutetium",     "Hafnium",     "Tantalum",    "Tungsten",      "Rhenium", 
  "Osmium",       "Iridium",     "Platinum",    "Gold",          "Mercury",
  "Thallium",     "Lead",        "Bismuth",     "Polonium",      "Astatine",
  "Radon",        "Francium",    "Radium",      "Actinium",      "Thorium", 
  "Protactinium", "Uranium",     "Neptunium",   "Plutonium",     "Americium", 
  "Curium",       "Berkelium",   "Californium", "Einsteinium",   "Fermium", 
  "Mendelevium",  "Nobelium",    "Lawrencium",  "Rutherfordium", "Dubnium",
  "Seaborgium",   "Bohrium",     "Hassium",     "Meitnerium",    "Darmstadtium",
  "Roentgenium"};

/* -------------------------------------------------------------------------- */

/* convert atomic number to atom name

   a - atomic number */
char *atom_name(unsigned a) {
  static char name[5];
  if (a && a<=ATOM_MAX_NUM)
    sprintf(name,"%s",atom_names[a-1]);
  else
    sprintf(name,"X");
  return(name);
  }

/* convert atomic number to atom full name

   a - atomic number */
char *atom_name_full(unsigned a) {
  static char name[30];
  if (a<=ATOM_MAX_NUM)
    sprintf(name,"%s",atom_full_names[a-1]);
  else
    sprintf(name,"Unknown");
  return(name);
  }

/* convert PDB atomic symbol to atom name

   s - atom name in PDB format */
char* atom_name_pdb(char *s) {
  static char name[5];
  int i = 0,j = 0,len;
  len = str_length(s);
  /* trim numbers */
  while (i<len && !isalpha(s[i]))
    i++;
  while (i<len && isalpha(s[i])) {
    name[j]=s[i];
    i++;
    j++;
    }
  name[j]='\0';
  /* special cases */
  if      (str_compare(name,"A"))   sprintf(name,"Au");
  else if (str_compare(name,"CA"))  sprintf(name,"C");
  else if (str_compare(name,"CB"))  sprintf(name,"C");
  else if (str_compare(name,"CC"))  sprintf(name,"C");
  else if (str_compare(name,"CD"))  sprintf(name,"C");
  else if (str_compare(name,"CE"))  sprintf(name,"C");
  else if (str_compare(name,"CG"))  sprintf(name,"C");
  else if (str_compare(name,"CH"))  sprintf(name,"C");
  else if (str_compare(name,"CM"))  sprintf(name,"C");
  else if (str_compare(name,"CZ"))  sprintf(name,"C");
  else if (str_compare(name,"CAA")) sprintf(name,"C");
  else if (str_compare(name,"CAB")) sprintf(name,"C");
  else if (str_compare(name,"CAC")) sprintf(name,"C");
  else if (str_compare(name,"CAD")) sprintf(name,"C");
  else if (str_compare(name,"CBA")) sprintf(name,"C");
  else if (str_compare(name,"CBB")) sprintf(name,"C");
  else if (str_compare(name,"CBC")) sprintf(name,"C");
  else if (str_compare(name,"CBD")) sprintf(name,"C");
  else if (str_compare(name,"CCB")) sprintf(name,"C");
  else if (str_compare(name,"CCB")) sprintf(name,"C");
  else if (str_compare(name,"CCD")) sprintf(name,"C");
  else if (str_compare(name,"CCF")) sprintf(name,"C");
  else if (str_compare(name,"CCH")) sprintf(name,"C");
  else if (str_compare(name,"CGA")) sprintf(name,"C");
  else if (str_compare(name,"CGD")) sprintf(name,"C");
  else if (str_compare(name,"CHA")) sprintf(name,"C");
  else if (str_compare(name,"CHB")) sprintf(name,"C");
  else if (str_compare(name,"CHC")) sprintf(name,"C");
  else if (str_compare(name,"CHD")) sprintf(name,"C");
  else if (str_compare(name,"CMA")) sprintf(name,"C");
  else if (str_compare(name,"CMB")) sprintf(name,"C");
  else if (str_compare(name,"CMC")) sprintf(name,"C");
  else if (str_compare(name,"CMD")) sprintf(name,"C");
  else if (str_compare(name,"Hb"))  sprintf(name,"H");
  else if (str_compare(name,"HA"))  sprintf(name,"H");
  else if (str_compare(name,"HB"))  sprintf(name,"H");
  else if (str_compare(name,"HC"))  sprintf(name,"H");
  else if (str_compare(name,"HD"))  sprintf(name,"H");
  else if (str_compare(name,"HE"))  sprintf(name,"H");
  else if (str_compare(name,"HG"))  sprintf(name,"H");
  else if (str_compare(name,"HH"))  sprintf(name,"H");
  else if (str_compare(name,"HM"))  sprintf(name,"H");
  else if (str_compare(name,"HN"))  sprintf(name,"H");
  else if (str_compare(name,"HO"))  sprintf(name,"H");
  else if (str_compare(name,"HP"))  sprintf(name,"H");
  else if (str_compare(name,"HS"))  sprintf(name,"H");
  else if (str_compare(name,"HT"))  sprintf(name,"H");
  else if (str_compare(name,"HV"))  sprintf(name,"H");
  else if (str_compare(name,"HW"))  sprintf(name,"H");
  else if (str_compare(name,"HZ"))  sprintf(name,"H");
  else if (str_compare(name,"HAA")) sprintf(name,"H");
  else if (str_compare(name,"HAB")) sprintf(name,"H");
  else if (str_compare(name,"HAC")) sprintf(name,"H");
  else if (str_compare(name,"HAD")) sprintf(name,"H");
  else if (str_compare(name,"HAM")) sprintf(name,"H");
  else if (str_compare(name,"HBA")) sprintf(name,"H");
  else if (str_compare(name,"HBB")) sprintf(name,"H");
  else if (str_compare(name,"HBC")) sprintf(name,"H");
  else if (str_compare(name,"HBD")) sprintf(name,"H");
  else if (str_compare(name,"HBM")) sprintf(name,"H");
  else if (str_compare(name,"HCB")) sprintf(name,"H");
  else if (str_compare(name,"HCF")) sprintf(name,"H");
  else if (str_compare(name,"HCM")) sprintf(name,"H");
  else if (str_compare(name,"HDM")) sprintf(name,"H");
  else if (str_compare(name,"HGM")) sprintf(name,"H");
  else if (str_compare(name,"HMA")) sprintf(name,"H");
  else if (str_compare(name,"HMB")) sprintf(name,"H");
  else if (str_compare(name,"HMC")) sprintf(name,"H");
  else if (str_compare(name,"HMD")) sprintf(name,"H");
  else if (str_compare(name,"HOA")) sprintf(name,"H");
  else if (str_compare(name,"HOB")) sprintf(name,"H");
  else if (str_compare(name,"HOD")) sprintf(name,"H");
  else if (str_compare(name,"HOF")) sprintf(name,"H");
  else if (str_compare(name,"HVC")) sprintf(name,"H");
  else if (str_compare(name,"HVT")) sprintf(name,"H");
  else if (str_compare(name,"NA"))  sprintf(name,"N");
  else if (str_compare(name,"NB"))  sprintf(name,"N");
  else if (str_compare(name,"NC"))  sprintf(name,"N");
  else if (str_compare(name,"ND"))  sprintf(name,"N");
  else if (str_compare(name,"NE"))  sprintf(name,"N");
  else if (str_compare(name,"NH"))  sprintf(name,"N");
  else if (str_compare(name,"NZ"))  sprintf(name,"N");
  else if (str_compare(name,"Ob"))  sprintf(name,"O");
  else if (str_compare(name,"OB"))  sprintf(name,"O");
  else if (str_compare(name,"OD"))  sprintf(name,"O");
  else if (str_compare(name,"OE"))  sprintf(name,"O");
  else if (str_compare(name,"OG"))  sprintf(name,"O");
  else if (str_compare(name,"OH"))  sprintf(name,"O");
  else if (str_compare(name,"OT"))  sprintf(name,"O");
  else if (str_compare(name,"OTi")) sprintf(name,"O");
  else if (str_compare(name,"OW"))  sprintf(name,"O");
  else if (str_compare(name,"OXT")) sprintf(name,"O");
  else if (str_compare(name,"SD"))  sprintf(name,"S");
  else if (str_compare(name,"SG"))  sprintf(name,"S");
  else if (str_compare(name,"Tib")) sprintf(name,"Ti");
  else if (str_compare(name,"TiO")) sprintf(name,"Ti");
  else if (str_compare(name,"TiS")) sprintf(name,"Ti");
  return(name);
  }

/* convert atomic name to PDB position style

   s - atom name in PDB format */
char* atom_name_to_pdb(char *s) {
  static char name[5];
  char *t;
  unsigned n;
  /* trim original name */
  t = str_copy_new(s);
  str_trim(t);
  n = str_length(t);
  /* set position */
  if (n>=4)
    sprintf(name,"%4.4s",t);
  else if (n>=3) {
    if ((t[0]>='A' && t[0]<='Z') || (t[0]>='a' && t[0]<='z'))
      sprintf(name," %s",t);
    else
      sprintf(name,"%s ",t);
    }
  else if (n>=2) {
    if ((t[0]>='A' && t[0]<='Z') || (t[0]>='a' && t[0]<='z'))
      sprintf(name," %s ",t);
    else
      sprintf(name,"%s  ",t);
    }
  else 
    sprintf(name," %s  ",t);
  return(name);
  }

/* -------------------------------------------------------------------------- */

/* convert atom name to atomic number

   nme - atom name */
unsigned atom_num(char *name) {
  unsigned i;
  char str[5];
  str_Lowcase_copy(name,str);
  for (i=0; i<ATOM_MAX_NUM; i++)
    if (str_compare(str,atom_names[i]))
      break;
  if (i<ATOM_MAX_NUM)
    return(i+1);
  return(0);
  }

/* convert atom name in PDB format to atomic number

   s - atom name */
unsigned atom_num_pdb(char *s) {
  return(atom_num(atom_name_pdb(s)));
  }

/* convert atomic mass to atomic number
 
   m - atomic mass */
unsigned atom_num_mass(double m) {
  unsigned i,n = 0;
  double v,d = 9.9e+99;
  for (i=0; i<ATOM_MAX_NUM; i++) {
    v = fabs(atom_mass(i+1)-m);
    if (v<d) {
      d = v;
      n = i+1;
      }
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
