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
#include "pw92eps.hpp"

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template <typename num> static num phi(const densvars<num> & d) {
  return pow(2.0, -1.0 / 3.0) * d.n_m13 * d.n_m13 * (sqrt(d.a_43) + sqrt(d.b_43));
}

template <typename num> static num energy(const densvars<num> & d) {
  using xcfun_constants::param_gamma;
  const parameter beta = 0.052;
  const parameter alpha = 1.0;
  num bg = beta / param_gamma;
  num eps = pw92eps::pw92eps(d);
  num u = phi(d);
  num u3 = pow3(u);
  // d2 is t^2
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * d.gnn /
           (u * u * pow(d.n, 7.0 / 3.0));
  num tt = pow(d2, 0.5); // this is t
  num v = tt * u * pow(d.r_s / 3.0, -1.0 / 6.0);
  num v3 = pow3(v);
  //
  // The term containing abs(zeta) gives problems with the
  // automatic derivatives. Here it is implemented substituting
  // abs(zeta)^omega with a polinomial fit. (see zvpbesolc.cpp
  // for more details).
  num zw = (0.462757 + 1.30129 * d.zeta * d.zeta -
            1.59546 * d.zeta * d.zeta * d.zeta * d.zeta +
            1.19635 * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta -
            0.36519 * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta *
                d.zeta) *
           d.zeta * d.zeta * d.zeta * d.zeta;
  //
  num ff = exp(-alpha * v3 * zw);
  num A = bg / expm1(-eps / (param_gamma * u3));
  num d2A = d2 * A;
  num H = param_gamma * u3 * log(1 + bg * d2 * (1 + d2A) / (1 + d2A * (1 + d2A)));
  return d.n * (eps + ff * H);
}

FUNCTIONAL(XC_ZVPBEINTC) = {"zvPBEint correlation Functional",
                            "zvPBEint correlation Functional\n"
                            "L.A. Constantin, E. Fabiano, F. Della Sala ,\n"
                            "      J. Chem. Phys. 137, 194105 (2012) .\n"
                            "Implemented by Eduardo Fabiano\n",
                            XC_DENSITY | XC_GRADIENT,
                            ENERGY_FUNCTION(energy)};
