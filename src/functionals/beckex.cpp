#include "functional.hpp"
#include "constants.hpp"

template<class num>
static num becke_alpha(const num &na, const num &gaa)
{
  const parameter c = pow(81/(4*M_PI),1.0/3.0)/2;
  const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num lda = -c*na43;
  num chi2 = gaa*pow(na,-8.0/3.0);
  num b88 = -(d*na43*chi2)/(1+6*d*sqrtx_asinh_sqrtx(chi2));
  return lda + b88;
}

// Lee, Lee and Paars kinetic energy functional.
// Phys. Rev. A 44, 768 (1991).
// new g = 0.003215546875
template<class num>
 static num llp91_k(const num &na, const num &gaa)
{
 using xc_constants::CF;
 const parameter c = pow(2.0,2.0/3.0)*CF;
 const parameter a = 4.4188E-3;
 const parameter g = 0.0253;
 num na53 = pow(na,5.0/3.0);
 num lda = c*na53;
 num chi2 = gaa*pow(na,-8.0/3.0);
 num gga = (a*c*na53*chi2)/(1+g*sqrtx_asinh_sqrtx(chi2));
 return lda + gga;
}

template<class num>
 static num llp91_ks(const num &na, const num &gaa)
{
  using xc_constants::CF;
  const parameter c = pow(2.0,2.0/3.0)*CF;
  const parameter a = 4.4188E-3;
  const parameter g = 0.03215546875;
  num na53 = pow(na,5.0/3.0);
  num lda = c*na53;
  num chi2 = gaa*pow(na,-8.0/3.0);
  num gga = (a*c*na53*chi2)/(1+g*sqrtx_asinh_sqrtx(chi2));
  return lda + gga;
}

template<class num>
static num becke_corr(const num &na, const num &gaa)
{
  const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num chi2 = gaa*pow(na,-8.0/3.0);
  return -(d*na43*chi2)/(1+6*d*sqrtx_asinh_sqrtx(chi2));
}

// Short range becke exchange as used in camb3lyp If mu=0 this reduces
// to the standard beckex for which must be used instead.
// FIXME: The erf + something is basically the erf minus
// its asymptotic expansion. This is horribel for numerics,
// will have to code a special function. 
// As coded here the code will fail if mu = 0, in which case the
// regular beckex should be used.

template<class num>
static num becke_sr(parameter mu, const num &na, const num &gaa)
{
  const parameter cparam = pow(81/(4*M_PI),1.0/3.0)/2;
  const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num chi2 = gaa*pow(na,-8.0/3.0);
  num K = 2*(cparam + (d*chi2)/(1+6*d*sqrtx_asinh_sqrtx(chi2)));
  num a = mu*sqrt(K)/(6*sqrt(M_PI)*pow(na,1.0/3.0));
  num b = expm1(-1/(4*a*a));
  num c = 2*a*a*b + 0.5;
  return -0.5*na43*K*(1-8.0/3.0*a*(sqrt(M_PI)*erf(1/(2*a))+2*a*(b-c)));
}

template<class num>
static num becke_cam(parameter alpha, parameter beta, parameter mu, const num &na, const num &gaa)
{
  const parameter cparam = pow(81/(4*M_PI),1.0/3.0)/2;
  const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num chi2 = gaa*pow(na,-8.0/3.0);
  num K = 2*(cparam + (d*chi2)/(1+6*d*sqrtx_asinh_sqrtx(chi2)));
  num a = mu*sqrt(K)/(6*sqrt(M_PI)*pow(na,1.0/3.0));
  num b = expm1(-1/(4*a*a));
  num c = 2*a*a*b + 0.5;
  return -0.5*na43*K*(1-alpha-beta*8.0/3.0*a*(sqrt(M_PI)*erf(1/(2*a))+2*a*(b-c)));
}

template<class num>
static num beckex(const densvars<num> &d)
{
  return becke_alpha(d.a,d.gaa) + becke_alpha(d.b,d.gbb);
}

template<class num>
static num beckexcorr(const densvars<num> &d)
{
  return becke_corr(d.a,d.gaa) + becke_corr(d.b,d.gbb);
}

template<class num>
static num beckesrx(const densvars<num> &d)
{
  parameter mu = d.get_param(XC_RANGESEP_MU);
  return becke_sr(mu,d.a,d.gaa) + becke_sr(mu,d.b,d.gbb);
}

template<class num>
static num beckecamx(const densvars<num> &d)
{
  parameter mu = d.get_param(XC_RANGESEP_MU);
  parameter alpha = d.get_param(XC_CAM_ALPHA);
  parameter beta = d.get_param(XC_CAM_BETA);
  return becke_cam(alpha,beta,mu,d.a,d.gaa) + becke_cam(alpha,beta,mu,d.b,d.gbb);
}

template<class num>
static num llp91k(const densvars<num> &d)
{
  return llp91_k(d.a, d.gaa) + llp91_k(d.b, d.gbb);
}

template<class num>
static num llp91ks(const densvars<num> &d)
{
  return llp91_ks(d.a, d.gaa) + llp91_ks(d.b, d.gbb);
}


FUNCTIONAL(XC_BECKEX) = {
  "Becke 88 exchange",
  "Becke 88 exchange including Slater part\n"
  "A.D. Becke, Density-functional exchange-energy approximation\n"
  "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
  "Implemented by Ulf Ekstrom\n"
  "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(beckex)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  {	-0.277987329958E+03,

	-0.385951846654E+01,
	-0.381309494319E+01,
	-0.172434478018E-04,
	0.000000000000E+00,
	-0.173712338362E-04,

	-0.441426807406E-01,
	0.000000000000E+00,
	0.201415922856E-06,
	0.000000000000E+00,
	0.000000000000E+00,
	-0.447245742260E-01,
	0.000000000000E+00,
	0.000000000000E+00,
	0.195961359539E-06,
	0.700742719647E-11,
	0.000000000000E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.718678968862E-11}
};
 
FUNCTIONAL(XC_BECKECORRX) = {
  "Becke 88 exchange correction",
  "Becke 88 exchange not including Slater part\n"
  "A.D. Becke, Density-functional exchange-energy approximation\n"
  "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
  "Implemented by Ulf Ekstrom\n"
  "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
  XC_DENSITY|XC_GRADIENT,
  ENERGY_FUNCTION(beckexcorr)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  { 0.39e+02,
    0.38e+02,
    0.81e+06,
    0.82e+06,
    0.82e+06},
  {
//     radovan: reference data obtained from *.c implementation in DIRAC
    -3.603918211981e+01, // 00000
    
    3.479609002901e-01, // 10000
    3.581112448092e-01, // 01000
    -1.724344780183e-05, // 00100
    0.000000000000e+00, // 00001
    -1.737123383621e-05, // 00010
    
    -8.181318630937e-03, // 20000
    0.000000000000e+00, // 11000
    2.014159228564e-07, // 10100
    0.000000000000e+00, // 10010
    0.000000000000e+00, // 10001
    -8.135046261131e-03, // 02000
    0.000000000000e+00, // 01100
    0.000000000000e+00, // 0100
    1.959613595393e-07, // 01010
    7.0074271964711398e-12, // radovan: i got this using xcfun
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    0.0000000000000000,
    7.1867896886212297e-12, // radovan: i got this using xcfun
  }
};

FUNCTIONAL(XC_LLP91K) = {
  "LLP91 kinetic energy functional",
  "LLP91 kinetic energy functional, Implemented by Moritz Bensberg\n"
  "Not an original part of xcFun.\n"
  "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
  XC_DENSITY|XC_GRADIENT,
  ENERGY_FUNCTION(llp91k)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  // Reference data obtained with working implementation in libXC
  // and compared to ADF results
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  { 4.584610965220e+03,

  8.418518384099e+01,
  8.268894675419e+01,
  2.798552163230e-04,
  0.000000000000E+00,
  2.794796623408e-04,

  1.486757132659e+00,
  0.000000000000E+00,
  -8.611445739344e-07,
  0.000000000000E+00,
  0.000000000000E+00,
  1.490091095937e+00,
  0.000000000000E+00,
  0.000000000000E+00,
  -6.848310104379e-07,
  -1.140141564535e-10,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  -1.159097944867e-10}
};

FUNCTIONAL(XC_BECKESRX) = {
  "Short range Becke 88 exchange",
  "Short range Becke 88 exchange, Implemented by Ulf Ekstrom\n"
  "Uses XC_RANGESEP_MU\n",
  XC_DENSITY|XC_GRADIENT,
  ENERGY_FUNCTION(beckesrx) };

FUNCTIONAL(XC_BECKECAMX) = {
  "CAM Becke 88 exchange",
  "CAM Becke 88 exchange, Implemented by Elisa Rebolini\n"
  "Uses XC_RANGESEP_MU\n",
  XC_DENSITY|XC_GRADIENT,
  ENERGY_FUNCTION(beckecamx) };

FUNCTIONAL(XC_LLP91KS) = {
  "LLP91 kinetic energy functional optimized for subsystem DFT",
  "LLP91 kinetic energy functional optimized for subsystem DFT\n"
  "Implemented by Moritz Bensberg\n"
  "Not an original part of xcFun.\n"
  "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
  XC_DENSITY|XC_GRADIENT,
  ENERGY_FUNCTION(llp91ks)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  // Reference data obtained with working implementation in xcFun
  // and compared to ADF results
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  { 4.538559033121e+03,

  8.510232713773e+01,
  8.365402093090e+01,
  2.460242017894e-04,
  0.000000000000E+00,
  2.446688531666e-04,

  1.448955945331e+00,
  0.000000000000E+00,
  1.044457826794e-07,
  0.000000000000E+00,
  0.000000000000E+00,
  1.451394306499e+00,
  0.000000000000E+00,
  0.000000000000E+00,
  2.818736159988e-07,
  -1.157859200527e-10,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  -1.167896572749e-10}
};

