// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MODIFIEDSPHERICALBESSEL1STKIND_H
#define INPSIGHTS_MODIFIEDSPHERICALBESSEL1STKIND_H

#include <vector>

namespace ModifiedSphericalBessel1stKind{

    std::vector<double> evaluateToMaxDegree(unsigned maxDegree, double x,
                                            double xZeroThreshold = std::numeric_limits<double>::epsilon(),
                                            double sphZeroThreshold = std::numeric_limits<double>::epsilon());

    std::vector<double> differentiateToMaxDegree(unsigned maxDegree, double x,
                                                 double xZeroThreshold = std::numeric_limits<double>::epsilon(),
                                                 double sphZeroThreshold = std::numeric_limits<double>::epsilon());
}

#endif //INPSIGHTS_MODIFIEDSPHERICALBESSEL1STKIND_H
