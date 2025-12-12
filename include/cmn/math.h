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

#ifndef ZF_LIB_CMN_MATH_H
#define ZF_LIB_CMN_MATH_H

/* -------------------------------------------------------------------------- */

#define MATH_3f2r7         1.13389341902768168164 /* 3/sqrt(7) */
#define MATH_3f2m2r2       1.06066017177982128660 /* 3/(2*sqrt(2)) */
#define MATH_3f2m2r7       0.56694670951384084082 /* 3/(2*sqrt(7)) */
#define MATH_3f2m2r14      0.40089186286863657702 /* 3/(2*sqrt(14)) */
#define MATH_5f2r21        1.09108945117996190633 /* 5/sqrt(21) */
#define MATH_5f2m2r2       1.76776695296636881100 /* 5/(2*sqrt(2)) */
#define MATH_2r2           1.41421356237309504880 /* sqrt(2) */
#define MATH_2r2f2         0.70710678118654752440 /* sqrt(2)/2 */
#define MATH_2r3           1.73205080756887729352 /* sqrt(3) */
#define MATH_2r3f2         0.86602540378443864676 /* sqrt(3)/2 */
#define MATH_2r3f2r2       1.22474487139158904909 /* sqrt(3)/sqrt(2) */
#define MATH_2r3f2m2r2     0.61237243569579452454 /* sqrt(3)/(2*sqrt(2)) */
#define MATH_2r3f2m2r5     0.38729833462074168851 /* sqrt(3)/(2*sqrt(5)) */
#define MATH_2r3f2m2r10    0.27386127875258305672 /* sqrt(3)/(2*sqrt(10)) */
#define MATH_2r5           2.23606797749978969640 /* sqrt(5) */
#define MATH_2r5f2         1.11803398874989484820 /* sqrt(5)/2 */
#define MATH_2r5f4         0.55901699437494742410 /* sqrt(5)/4 */
#define MATH_2r5f8         0.27950849718747371205 /* sqrt(5)/8 */
#define MATH_2r5f2r3       1.29099444873580562839 /* sqrt(5)/sqrt(3) */
#define MATH_2r5f2r6       0.91287092917527685576 /* sqrt(5)/sqrt(6) */
#define MATH_2r5f2m2r2     0.79056941504209483299 /* sqrt(5)/(2*sqrt(2)) */
#define MATH_2r5f2m2r3     0.64549722436790281419 /* sqrt(5)/(2*sqrt(3)) */
#define MATH_2r5f2m2r7     0.42257712736425828875 /* sqrt(5)/(2*sqrt(7)) */
#define MATH_2r5f4m2r6     0.22821773229381921394 /* sqrt(5)/(4*sqrt(6)) */
#define MATH_2r5f4m2r7     0.21128856368212914437 /* sqrt(5)/(4*sqrt(7)) */
#define MATH_2r5f8m2r3     0.16137430609197570354 /* sqrt(5)/(8*sqrt(3)) */
#define MATH_2r6f2r5       1.09544511501033222691 /* sqrt(6)/sqrt(5) */
#define MATH_2r10f2r7      1.19522860933439363996 /* sqrt(10)/sqrt(7) */
#define MATH_2r15f8        0.48412291827592711064 /* sqrt(15)/8 */
#define MATH_2r15f4m2r7    0.36596252735569994225 /* sqrt(15)/(4*sqrt(7)) */
#define MATH_2r35f4m2r3    0.85391256382996653193 /* sqrt(35)/(4*sqrt(3)) */
#define MATH_2r35f8        0.73950997288745200532 /* sqrt(35)/8 */
#define MATH_2r35f8m2r2    0.52291251658379721748 /* sqrt(35)/(8*sqrt(2)) */
#define MATH_3m2r3f4       1.29903810567665797014 /* 3*sqrt(3)/4 */
#define MATH_3m2r3f2r35    0.87831006565367986141 /* 3*sqrt(3)/sqrt(35) */
#define MATH_3m2r3f2m2r7   0.98198050606196571569 /* 3*sqrt(3)/(2*sqrt(7)) */
#define MATH_3m2r3f4m2r35  0.21957751641341996535 /* 3*sqrt(3)/(4*sqrt(35)) */
#define MATH_5m2r3f4m2r2   1.53093108923948631137 /* 5*sqrt(3)/(4*sqrt(2)) */
#define MATH_3m2r5f2m2r7   1.26773138209277486626 /* 3*sqrt(5)/(2*sqrt(7)) */
#define MATH_3m2r5f2m2r14  0.89642145700079522997 /* 3*sqrt(5)/(2*sqrt(14)) */
#define MATH_3m2r7f8m2r2   0.70156076002011400979 /* 3*sqrt(7)/(8*sqrt(2)) */
#define MATH_5m2r7f8m2r2   1.16926793336685668299 /* 5*sqrt(7)/(8*sqrt(2)) */

#define MATH_PI            3.14159265358979323844 /* Pi */
#define MATH_2PI           6.28318530717958647688 /* 2*Pi */
#define MATH_3PI           9.42477796076937971532 /* 3*Pi */
#define MATH_PIf2          1.57079632679489661922 /* Pi/2 */
#define MATH_PIf3          1.04719755119659774614 /* Pi/3 */
#define MATH_PIs2          9.86960440108935861869 /* Pi^2 */
#define MATH_PIs3         31.00627668029982017480 /* Pi^3 */
#define MATH_2r8PIs3      15.74960994572241974429 /* sqrt(8*Pi^3) */

#define MATH_1fPI          0.31830988618379067154 /* 1/Pi */
#define MATH_2fPI          0.63661977236758134308 /* 2/Pi */
#define MATH_3fPI          0.95492965855137201462 /* 3/Pi */
#define MATH_2r1fPI        0.56418958354775628695 /* sqrt(1/Pi) */
#define MATH_2r2fPI        0.79788456080286535588 /* sqrt(2/Pi) */
#define MATH_2r3fPI        0.97720502380583984317 /* sqrt(3/Pi) */
#define MATH_3r1fPI        0.68278406325529568147 /* (1/Pi)^(1/3) */
#define MATH_3r2fPI        0.86025401382809962534 /* (2/Pi)^(1/3) */
#define MATH_3r3fPI        0.98474502184269654118 /* (3/Pi)^(1/3) */

#define MATH_INV_03        0.33333333333333333333  /* 1/3 */
#define MATH_INV_06        0.16666666666666666667  /* 1/6 */
#define MATH_INV_07        0.14285714285714285714  /* 1/7 */
#define MATH_INV_09        0.11111111111111111111  /* 1/9 */
#define MATH_INV_11        0.09090909090909090909  /* 1/11 */
#define MATH_INV_12        0.08333333333333333333  /* 1/12 */
#define MATH_INV_13        0.07692307692307692308  /* 1/13 */
#define MATH_INV_14        0.07142857142857142857  /* 1/14 */
#define MATH_INV_15        0.06666666666666666667  /* 1/15 */
#define MATH_INV_17        0.05882352941176470588  /* 1/17 */
#define MATH_INV_18        0.05555555555555555556  /* 1/18 */
#define MATH_INV_19        0.05263157894736842105  /* 1/19 */
#define MATH_INV_21        0.04761904761904761904  /* 1/21 */
#define MATH_INV_22        0.04545454545454545454  /* 1/22 */
#define MATH_INV_23        0.04347826086956521739  /* 1/23 */
#define MATH_INV_24        0.04166666666666666666  /* 1/24 */

/* -------------------------------------------------------------------------- */

/* return greater of the two values (int) */
int math_imax(int, int);
/* return smaller of the two values (int)*/
int math_imin(int, int);
/* return greater of the two values (unsigned) */
unsigned math_umax(unsigned, unsigned);
/* return smaller of the two values (unsigned)*/
unsigned math_umin(unsigned, unsigned);
/* return greater of the two values (double) */
double math_fmax(double, double);
/* return smaller of the two values (double) */
double math_fmin(double, double);

/* swap values between two variables (integers) */
void math_iswap(int*, int*);
/* swap values between two variables (unsigned int-s) */
void math_uswap(unsigned*, unsigned*);
/* swap values between two variables (short int-s) */
void math_siswap(short*, short*);
/* swap values between two variables (doubles) */
void math_fswap(double*, double*);

/* calculate length of hypotenuse by pythagorean theorem */
double math_fpyth(double, double);

/* evaluate factorial function - singed integer version */
double math_fact_if(int n);
/* evaluate factorial function */
double math_fact_f(unsigned);
/* evaluate logarithm of factorial function */
double math_fact_ln_f(unsigned);
/* evaluate odd factorial function */
double math_fact_odd_f(unsigned);

/* evaluate binomical function (real version) */
double math_binom_f(unsigned, unsigned);
/* evaluate binomical function (integer version) */
unsigned math_binom_u(unsigned, unsigned);

/* calculate arc cosine (argument truncation to [-1,+1] interval) */
double math_acos(double c);

/* evaluate error function */
double math_erf(double);
/* evaluate complementary error function */
double math_erfc(double);
/* evaluate inversion of complementary error function */
double math_erfc_inv(double);
/* evaluate complementary error function by Chebyshev interpolation */
double math_erfc_chebyshev(double);

/* evaluate logarithm of gamma function */
double math_gamma_ln(double);
/* evaluate incomplete gamma function p(a,x) */
double math_igamma_p(double, double);
/* evaluate inverse incomplete gamma function p[-1](a,x) */
double math_igamma_p_inv(double, double);
/* evaluate incomplete gamma function q(a,x) */
double math_igamma_q(double, double);
/* evaluate beta function */
double math_beta(double, double);

/* evaluate integral Boys function (core of GTO integrals) */
double math_boys(unsigned, double);

/* calculate arithmetic average and variance from array of values */
double math_avrg(double*, long unsigned, double*);

/* evaluate normalized gaussian function */
double math_gauss_fce(double, double, double);
/* fit values to normalized gaussian function */
void math_gauss_fit(double*, double*, unsigned, double*, double*);

/* evaluate normalized cauchy function */
double math_cauchy_fce(double, double, double);

/* calculate histogram from given values */
unsigned long* math_hist_calc(double*, long unsigned, double, double,
  unsigned, double*);
/* create normalized histogram from number of counts */
double *math_hist_norm(unsigned long*, unsigned, double);

/* converts spherical coordinates (r,p,t) to cartesian (x,y,z) */
void math_rpt_xyz(double*, double, double, double);
/* calculate lenght of triangle sides on the unit sphere */
void math_rpt_area_3p_sides(double*, double*, double*,
  double*, double*, double*);
/* calculate angles between triangle sides on the unit sphere */
void math_rpt_area_3p_angles(double, double, double, 
  double*, double*, double*);
/* calculate are selected by 3 points on the unit sphere */
double math_rpt_area_3p(double*, double*, double*);
/* calculate coordinates of spherical triangle centroid */
void math_rpt_area_3p_center(double*, double*, double*, double*);

/* integrate function by trapezoidal rule */
double math_integ_trapz(double*, double, double, unsigned);
/* integrate function by trapezoidal rule - 1/N3 accuracy */
double math_integ_trapz3n(double*, double, double, unsigned);
/* integrate function by trapezoidal rule - 1/N4 accuracy */
double math_integ_trapz4n(double*, double, double, unsigned);
/* integrate function by Simpson's rule */
double math_integ_simpson(double*, double, double, unsigned);
/* calculate spherical integral from given function 
   (longitude/latitude grid with centroid rule) */
double math_integ_sphere_llc (double, double,
  double f(double*, void*), void*);

/* Lowdin orthogonalization of m n-dimensional vectors */
void math_orth_lowdin_norm(double**, unsigned, unsigned);
/* Gramm-Schmidth orthogonalization of m n-dimensional vectors */
void math_orth_gs(double**, unsigned, unsigned);
/* Gramm-Schmidth orthonormalization of m n-dimensional vectors */
void math_orth_gs_norm(double**, unsigned, unsigned);
/* check orthogonality of given vectors */
short math_orth_check(double**, unsigned, unsigned);
/* check orthonormality of given vectors */
short math_orth_check_norm(double**, unsigned, unsigned);
/* print table of vector products and return if the vectors are orthogonal */
short math_orth_print(double**, unsigned, unsigned);
/* print table of vector products and return if the vectors are orthonormal */
short math_orth_print_norm(double**, unsigned, unsigned);

/* solve set of linear equations by Gauss elimination */
void math_eqset_gauss(double**, double*, unsigned);

/* find interval to which belongs the interpolation point */
short math_ipol_locate(double, double*, long unsigned,
  long unsigned*);
/* linear interpolation between two points */
double math_ipol_linear(double, double, double, double, double);
/* polynomial interpolation between two points */
double math_ipol_polynomial(double, double*, double*,
  long unsigned);
/* cubic-spline interpolation between two points */
double math_ipol_spline(double, double*, double*,
  long unsigned, long unsigned);

/* -------------------------------------------------------------------------- */

#endif
