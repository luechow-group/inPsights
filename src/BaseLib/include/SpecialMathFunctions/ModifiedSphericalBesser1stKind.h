// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H
#define INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H

#include <vector>

struct ModifiedSphericalBessel1stKind
{
    ModifiedSphericalBessel1stKind(int degree);
    void evaluate(double r, bool differentiate);

    static std::vector<double> eval(int degree, double r);
    static constexpr double RADZERO = 1e-10;
    static constexpr double SPHZERO = 1e-4;

    int _degree;
    std::vector<double> _in;
    std::vector<double> _din;
};

#endif //INPSIGHTS_MODIFIEDSPHERICALBESSER1STKIND_H
