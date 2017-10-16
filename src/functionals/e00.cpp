#include "functional.hpp"
#include "pw9xx.hpp"

//Kinetic energy functional proposed by M. Ernzerhof in 2000
//THEOCHEM 501-502 (2000) 59-64
template<class num>
static num energy_e00(const num &na, const num &gaa)
{
  using pw91_like_x_internal::S2;
  using xc_constants::CF;
  num st2 = S2(na, gaa);
  num na53 = pow(na,5.0/3.0);
  num lda = CF*na53;
  num e00 = lda * ((135+28*st2+5*st2*st2)/(135+3*st2));
  return e00;
}

template <class num>
static num energyE00(const densvars<num> &d)
{
  return energy_e00(d.a,d.gaa) + energy_e00(d.b,d.gbb);
}

FUNCTIONAL(XC_E00) = {
  "M. Ernzerhofs Kinetic Energy Functional",
  "M. Ernzerhofs Kinetic Energy Functional\n"
  "M. Ernzerhof, THEOCHEM 501-502 (2000) 59-64\n"
  "Implemented by Moritz Bensberg\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energyE00)
};



