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
#include "cmn/math.h"
#include "cmn/message.h"

#define GM_EPS     1.0E-15  /* gamma function convergence criterion */
#define GM_FPMIN   1.0E-50  /* gamma function floating point minimum */
#define GM_INV_EPS 1.0E-08  /* inverse gamma function convergence criterion */

/* -------------------------------------------------------------------------- */

/* internal data for gamma function */
static unsigned gm_ngau=18;
static double gm_y[18]={
  0.0021695375159141994, 0.011413521097787704,
  0.0279723089503021160, 0.051727015600492421,
  0.0825022254843409410, 0.120070199109602930,
  0.1641528330075247000, 0.214423769867793550,
  0.2705108284064433600, 0.331998763414478870,
  0.3984323418640194300, 0.469319714073754830,
  0.5441360555665797300, 0.622327452880310770,
  0.7033150046559717400, 0.786499107683134470,
  0.8712638961906151700, 0.956981801526291420};
static double gm_w[18]={
  0.0055657196642445571, 0.012915947284065419,
  0.0201815152977353820, 0.027298621498568734,
  0.0342138107702995370, 0.040875750923643261,
  0.0472350834902655820, 0.053244713977759692,
  0.0588601442453247980, 0.064039797355015485,
  0.0687453238357364080, 0.072941885005653087,
  0.0765984106458706400, 0.079687828912071670,
  0.0821872667043397060, 0.084078218979661945,
  0.0853466857393387210, 0.085983275670394821};
static double gm_gln;

/* -------------------------------------------------------------------------- */

/* evaluate logarithm of gamma function */
double math_gamma_ln(double a) {
  double x,y,t,s;
  unsigned i;
  static double coeff[14] = {
    57.156235665862923500e+0,-59.597960355475491200e+0,
    14.136097974741747100e+0, -0.491913816097620199e+0,
     0.339946499848118887e-4,  0.465236289270485756e-4,
    -0.983744753048795646e-4,  0.158088703224912494e-3,
    -0.210264441724104883e-3,  0.217439618115212643e-3,
    -0.164318106536763890e-3,  0.844182239838527433e-4,
    -0.261908384015814087e-4,  0.368991826595316234e-5};
  if (a<=0.0)
    msg_error("invalid argument of incomplete gamma function",1);
  y = x = a;
  t = x+5.24218750000000000;
  t = (x+0.5)*log(t)-t;
  s = 0.999999999999997092;
  for (i=0; i<14; i++)
    s += coeff[i]/++y;
  return(t+log(2.5066282746310005*s/x));
  }

/* evaluate incomplete gamma function p(a,x) by series representation */
double math_igamma_p_series(double a, double x) {
  double ap = a,del,sum;
  gm_gln = math_gamma_ln(a);
  del = sum = 1.0/a;
  for (;;) {
    ap++;
    del *= x/ap;
    sum += del;
    if (fabs(del)<fabs(sum)*GM_EPS)
      return(sum*exp(-x+a*log(x)-gm_gln));
    }
  return(0.0);
  }

/* evaluate incomplete gamma function q(a,x) by continued fractions */
double math_igamma_q_cfrac(double a, double x) {
  double an,b,c,d,h,del;
  long int i;
  gm_gln = math_gamma_ln(a);
  b = x+1.0-a;
  c = 1.0/GM_FPMIN;
  d = 1.0/b;
  h = d;
  for (i=1;;i++) {
    an = (-i*(i-a));
    b += 2.0;
    d = (an*d+b);
    if (fabs(d)<GM_FPMIN)
      d = GM_FPMIN;
    c = (b+an/c);
    if (fabs(c)<GM_FPMIN)
      c = GM_FPMIN;
    d = 1.0/d;
    del = (d*c);
    h *= del;
    if (fabs(del-1.0)<=GM_EPS)
      break;
    }
  return(exp(-x+a*log(x)-gm_gln)*h);
  }

/* evaluate incomplete gamma function p(a,x) or q(a,x) by quadrature */
double math_igamma_p_approx(double a, double x, int p) {
  double t,xu,a1,lna1,sqrta1,sum,ans;
  unsigned i;
  a1 = a-1.0;
  lna1 = log(a1);
  sqrta1 = sqrt(a1);
  gm_gln = math_gamma_ln(a);
  if (x>a1)
    xu = math_fmax(a1+11.5*sqrta1,x+6.0*sqrta1);
  else
    xu = math_fmax(0.0,math_fmin(a1-7.5*sqrta1,x-5.0*sqrta1));
  sum=0;
  for (i=0; i<gm_ngau; i++) {
    t = x+(xu-x)*gm_y[i];
    sum += (gm_w[i]*exp(-(t-a1)+a1*(log(t)-lna1)));
    }
  ans = (sum*(xu-x)*exp(a1*(lna1-1.0)-gm_gln));
  return(p ? (ans>0.0 ? 1.0-ans : -ans) : (ans>-0.0 ? ans : 1.0+ans));
  }

/* evaluate incomplete gamma function p(a,x) */
double math_igamma_p(double a, double x) {
  if (x<0.0 || a<=0.0)
    msg_error("invalid arguments in incomplete gamma function p(a,x)",1);
  if (x==0.0)
    return(0.0);
  else if ((int)a>=100)
    return(math_igamma_p_approx(a,x,1));
  else if (x<(a+1.0))
    return(math_igamma_p_series(a,x));
  return(1.0-math_igamma_q_cfrac(a,x));
  }

/* evaluate inverse incomplete gamma function p[-1](a,x) */
double math_igamma_p_inv(double p, double a) {
  double a1,err,lna1,afac,pp,u,t,x;
  unsigned i;
  a1 = a-1;
  if (a<=0.0)
    msg_error("negative order of inverser incomplete gamma function p(p,a)",1);
  if (p>=1.0)
    return(math_fmax(100.0,a+100.0*sqrt(a)));
  if (p<=0.0)
    return(0.0);
  afac = 1.0;
  lna1 = 0.0;
  if (a>1.0) {
    lna1 = log(a1);
    afac = exp(a1*(lna1-1.0)-gm_gln);
    pp = (p<0.5 ? p : 1.0-p);
    t = sqrt(-2.0*log(pp));
    x = (2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481))-t;
    if (p<0.5)
      x = -x;
    x = math_fmax(1.0E-3,a*pow(1.0-1.0/(9.0*a)-x/(3.0*sqrt(a)),3));
    }
  else {
    t = 1.0-a*(0.253+a*0.12);
    if (p<t)
      x = pow(p/t,1.0/a);
    else
      x = 1.0-log(1.0-(p-t)/(1.0-t));
    }
  for (i=0; i<12; i++) {
    if (x<=0.0)
      return(0.0);
    err = math_igamma_p(a,x)-p;
    if (a>1.0) 
      t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
    else
      t = exp(-x+a1*log(x)-gm_gln);
    u = err/t;
    x -= (t=u/(1.0-0.5*math_fmin(1.0,u*((a-1.0)/x-1))));
    if (x<=0.0)
      x = 0.5*(x+t);
    if (fabs(t)<GM_INV_EPS*x)
      break;
  }
  return(x);
  }

/* evaluate incomplete gamma function q(a,x) */
double math_igamma_q(double a, double x) {
  if (x<0.0 || a<=0.0)
    msg_error("invalid arguments in incomplete gamma function q(a,x)",1);
  if (x==0.0)
    return(1.0);
  else if ((int)a>=100)
    return(math_igamma_p_approx(a,x,0));
  else if (x<(a+1.0))
    return(1.0-math_igamma_p_series(a,x));
  return(math_igamma_q_cfrac(a,x));
  }

/* -------------------------------------------------------------------------- */

/* evaluate beta function */
double math_beta(double z, double w) {
  return(exp(math_gamma_ln(z)+math_gamma_ln(w)-math_gamma_ln(z+w)));
  }

/* -------------------------------------------------------------------------- */
