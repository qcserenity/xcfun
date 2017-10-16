#include "functional.hpp"
#include "pw9xx.hpp"

// PBE2/PBE3/PBE4 kinetic energy Functionals by V.V. Karasiev et. al.
// multiplied by 2^(2/3) to fit the results obtained with ADF
// V.V. Karasiev, S.B. Trickey and Frank E. Harris, J. Comput.-Aid. Mat. Des. 13, 111 (2006)
template<class num>
static num energy_pbe2(const num &na, const num &gaa)
{
  using pw91_like_x_internal::S2;
  using xc_constants::CF;
  const parameter c1 =2.0309;
  const parameter a1 =0.2942;
  num st2 = S2(na, gaa);
  num na53 = pow(na,5.0/3.0);
  num lda = pow(2.0,2.0/3.0)*CF*na53;
  num t1 = c1*st2;
  num t2 = 1+a1*st2;
  num pbe2 = lda * (1+t1/t2);
  return pbe2;
}

template<class num>
static num energy_pbe3(const num &na, const num &gaa)
{
  using pw91_like_x_internal::S2;
  using xc_constants::CF;
  const parameter c1 =-3.7425;
  const parameter c2 =50.258;
  const parameter a1 =4.1355;
  num st2 = S2(na, gaa);
  num na53 = pow(na,5.0/3.0);
  num lda = pow(2.0,2.0/3.0)*CF*na53;
  num a1st2_2=(1+a1*st2)*(1+a1*st2);
  num st4 = st2*st2;
  num pbe3 = lda * (1+c1*st2/(1+a1*st2)+c2*st4/a1st2_2);
  return pbe3;
}

template<class num>
static num energy_pbe4(const num &na, const num &gaa)
{   
  using pw91_like_x_internal::S2;
  using xc_constants::CF;
  const parameter c1 =-7.2333;
  const parameter c2 =61.645;
  const parameter c3 =-93.683;
  const parameter a1 =1.7107;
  num st2 = S2(na, gaa);
  num na53 = pow(na,5.0/3.0);
  num lda = pow(2.0,2.0/3.0)*CF*na53;
  num a1st2_2=(1+a1*st2)*(1+a1*st2);
  num a1st2_3=a1st2_2*(1+a1*st2);
  num st4 = st2*st2;
  num st6 = st2*st2*st2;
  num pbe4 = lda * (1+c1*st2/(1+a1*st2)+c2*st4/a1st2_2+c3*st6/a1st2_3);
  return pbe4;
}  

//PBE2 Functional with optimized parameters for subsystem DFT
//Parameters optimized by Danny Schlüns
template<class num>
static num energy_pbe2S(const num &na, const num &gaa)
{
  using pw91_like_x_internal::S2;
  using xc_constants::CF;
  const parameter c1 =-7.9691;
  const parameter a1 =31.975108203125;
  num st2 = S2(na, gaa);
  num na53 = pow(na,5.0/3.0);
  num lda = pow(2.0,2.0/3.0)*CF*na53;
  num t1 = c1*st2;
  num t2 = 1+a1*st2;
  num pbe2 = lda * (1+t1/t2);
  return pbe2;
}

template <class num>
static num energyPBE2(const densvars<num> &d)
{
  return energy_pbe2(d.a,d.gaa) + energy_pbe2(d.b,d.gbb);
}

template <class num>
static num energyPBE3(const densvars<num> &d)
{
  return energy_pbe3(d.a,d.gaa) + energy_pbe3(d.b,d.gbb);
}

template <class num>
static num energyPBE4(const densvars<num> &d)
{
  return energy_pbe4(d.a,d.gaa) + energy_pbe4(d.b,d.gbb);
}

template<class num>
static num energyPBE2S(const densvars<num> &d)
{  
   return energy_pbe2S(d.a, d.gaa) + energy_pbe2S(d.b, d.gbb);
}  

FUNCTIONAL(XC_PBE2) = {
  "PBE2 Kinetic Energy Functional",
  "PBE2 Kinetic Energy Functional\n"
  "V.V. Karasiev, S.B. Trickey and Frank E. Harris, J. Comput.-Aid. Mat. Des. 13, 111 (2006)\n"
  "Implemented by Moritz Bensberg\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energyPBE2)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  // Reference data obtained with working implementation in xcFun
  // and compared to ADF results
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  { 9.229855135481e+03,

  5.345663118734e+01,
  5.258357928260e+01,
  2.606126017070e-03,
  0.000000000000E+00,
  2.593387987753e-03,

  1.000126954100e+00,
  0.000000000000E+00,
  -1.558887539356e-06,
  0.000000000000E+00,
  0.000000000000E+00,
  7.520330815643e-01,
  0.000000000000E+00,
  0.000000000000E+00,
  2.962711473607e-06,
  -1.178393242146e-09,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  -1.237486748666e-09}
};

FUNCTIONAL(XC_PBE3) = {
  "PBE3 Kinetic Energy Functional",
  "PBE3 Kinetic Energy Functional\n"
  "V.V. Karasiev, S.B. Trickey and Frank E. Harris, J. Comput.-Aid. Mat. Des. 13, 111 (2006)\n"
  "Implemented by Moritz Bensberg\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energyPBE3)
};

FUNCTIONAL(XC_PBE4) = {
  "PBE4 Kinetic Energy Functional",
  "PBE4 Kinetic Energy Functional\n"
  "V.V. Karasiev, S.B. Trickey and Frank E. Harris, J. Comput.-Aid. Mat. Des. 13, 111 (2006)\n"
  "Implemented by Moritz Bensberg\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energyPBE4)
};

FUNCTIONAL(XC_PBE2S) = {
  "PBE2 Kinetic Energy Functional opt. for Subsystem DFT",
  "PBE2 Kinetic Energy Correction opt. for Subsystem DFT\n"
  "V.V. Karasiev, S.B. Trickey and Frank E. Harris, J. Comput.-Aid. Mat. Des. 13, 111 (2006)||TODO\n"
  "Implemented by Moritz Bensberg\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energyPBE2S)
    XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  // Reference data obtained with working implementation in xcFun
  // and compared to ADF results
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  { 3.042343568778e+03,

  6.776522476991e+01,
  6.644407128122e+01,
  -2.381993907754e-05,
  0.000000000000E+00,
  -2.089700838806e-05,

  1.297851503863e+00,
  0.000000000000E+00,
  -2.518237100635e-06,
  0.000000000000E+00,
  0.000000000000E+00,
  1.296660305162e+00,
  0.000000000000E+00,
  0.000000000000E+00,
  -2.276087229565e-06,
  5.649591944551e-11,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  0.000000000000E+00,
  4.911051361808e-11}
};
