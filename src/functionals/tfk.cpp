/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
 *
 * This file is part of XCFun.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * XCFun library, see: <https://xcfun.readthedocs.io/>
 */

#include "constants.hpp"
#include "functional.hpp"

//  Thomas-Fermi kinetic energy functional

template<class num>
static num tfk_energy(const num &na)
{
  using xcfun_constants::CF;

  return CF*pow(2.0,2.0/3.0)*pow(na, 5.0/3.0);
}

template<class num>
static num tfk(const densvars<num> &d)
{
    return tfk_energy(d.a)+tfk_energy(d.b);
}

FUNCTIONAL(XC_TFK) = {"Thomas-Fermi Kinetic Energy Functional",
                      "Thomas-Fermi Kinetic Energy Functional\n"
                      "Implemented by Andre Gomes.\n",
                      XC_DENSITY,
                      ENERGY_FUNCTION(tfk) XC_A_B,
                      XC_PARTIAL_DERIVATIVES,
                      1,
                      1e-5,
                      {1., .8},
                      {7.64755771625168,
                       7.08107195949229,
                       7.08107195949229,
                       2.62261924425641,
                       2.62261924425641,
                       2.62261924425641}};
