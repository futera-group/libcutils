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

#ifndef ZF_LIB_CMN_UNITS_H
#define ZF_LIB_CMN_UNITS_H

/* -------------------------------------------------------------------------- */

/* Atomic units */
#define UNIT_AU_KB         3.166829638E-06  /* Boltzmann constant */
#define UNIT_AU_H          1.000000000E+00  /* Hartree */
#define UNIT_AU_EV         3.674930883E-02  /* eV */
#define UNIT_AU_KCAL       1.593601219E-03  /* kcal/mol */
#define UNIT_AU_KJM        3.808798326E-04  /* kJ/mol */
#define UNIT_AU_ME         1.000000000E+00  /* mass of electron */
#define UNIT_AU_AMU        1.822888538E+03  /* atomic mass unit */
#define UNIT_AU_B          1.000000000E+00  /* Bohr */
#define UNIT_AU_ANG        1.889726125E+00  /* Angstrom */
#define UNIT_AU_NM         1.889726125E+01  /* nanometer */
#define UNIT_AU_FS         4.134137334E+01  /* femtosecond */
#define UNIT_AU_C          1.370359996E+02  /* speed of light */
#define UNIT_AU_D          3.934303000E-01  /* Debye */
#define UNIT_AU_PLK        6.283185307E+00  /* Planck constant */
#define UNIT_AU_NA         6.022136700E+23  /* Avogadro constant */
#define UNIT_AU_AMP        1.509748928E+02  /* Ampere */

/* SI units */
#define UNIT_SI_KB         1.380648813E-23  /* Boltzmann constant */
#define UNIT_SI_H          4.359744177E-18  /* Hartree [J] */
#define UNIT_SI_EV         1.602176462E-19  /* eV */
#define UNIT_SI_KCAL       6.947693637E-21  /* kcal/mol */
#define UNIT_SI_KJM        1.660538632E-21  /* kJ/mol */
#define UNIT_SI_ME         9.109381880E-31  /* mass of electron [kg] */
#define UNIT_SI_AMU        1.660538782E-27  /* atomic mass unit [kg] */
#define UNIT_SI_B          5.291772108E-11  /* Bohr [m] */
#define UNIT_SI_ANG        1.000000000E-10  /* Angstrom [m] */
#define UNIT_SI_FS         1.000000000E-15  /* femtosecond */
#define UNIT_SI_C          2.997924580E+08  /* speed of light [m/s] */
#define UNIT_SI_PLK        6.626068963E-34  /* Planck constant */
#define UNIT_SI_NA         6.022136700E+23  /* Avogadro constant [mol] */

/* Conversion */

/* Length */

#define CONV_B_ANG         5.291772108E-01  /* Bohr -> Angstrom */

/* Energy */
#define CONV_H_EV          2.721139613E+01  /* H -> eV */
#define CONV_H_KCAL        6.275095600E+02  /* H -> kcal/mol */
#define CONV_H_KJM         2.625499999E+03  /* H -> kJ/mol */
#define CONV_CAL_J         4.184000000E+00  /* calorie -> Joule */
#define CONV_ERG_J         1.000000000E-07  /* erg -> Joule */

/* Mass */
#define CONV_ME_AMU        5.485798934E-04  /* electron mass -> atom mass */

/* Force */
#define CONV_Dyne_JM       1.000000000E-05  /* Dyne -> Joule metre */

/* Force constant */
#define CONV_mDyneA_HB2    0.064229619E+00  /* mDyne/A -> H/B2 */
#define CONV_mDyneA_KcalA2 1.439306972E+02  /* mDyne/A -> kcal/mol/A2 */

/* Dipole moment */
#define CONV_D_eA          2.081917990E-01  /* Debye -> eA */

/* Pressure */
#define CONV_ATM_Pa        1.013250000E+05  /* atm -> Pa */

/* -------------------------------------------------------------------------- */

#endif
