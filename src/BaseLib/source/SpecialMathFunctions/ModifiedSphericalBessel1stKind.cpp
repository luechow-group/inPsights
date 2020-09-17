// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include "SpecialMathFunctions/ModifiedSphericalBessel1stKind.h"


std::vector<double> ModifiedSphericalBessel1stKind::evaluateToMaxDegree(unsigned maxDegree, double x,
                                                                        double xZeroThreshold,
                                                                        double sphZeroThreshold) {
    std::vector<double> in(maxDegree + 1);

    if (x < xZeroThreshold) {
        in[0] = (1.0);
        std::fill(std::next(std::begin(in)), std::end(in), 0);
    } else {
        in[0] = (std::sinh(x) / x);

        if (maxDegree >= 1) {
            in[1] = (std::cosh(x) / x - std::sinh(x) / (x * x));

            for (unsigned n = 2; n <= maxDegree; ++n) {
                if (in[n - 1] < sphZeroThreshold)
                    in[n] = 0.0;
                else
                    in[n] = (in[n - 2] - (2 * (n - 1) + 1) / x * in[n - 1]);
            }
        }
    }
    return in;
}

std::vector<double> ModifiedSphericalBessel1stKind::differentiateToMaxDegree(unsigned int maxDegree, double x,
                                                                             double xZeroThreshold,
                                                                             double sphZeroThreshold) {
    auto in = ModifiedSphericalBessel1stKind::evaluateToMaxDegree(maxDegree, x, xZeroThreshold, sphZeroThreshold);

    std::vector<double> din(maxDegree + 1, 0.0);
    if (x < xZeroThreshold) {
        din[1] = 1. / 3.;
    } else {
        din[0] = in[1];
        for (unsigned n = 1; n <= maxDegree; ++n)
            din[n] = in[n - 1] - (n + 1.0) / x * in[n];
    }
    return din;
}