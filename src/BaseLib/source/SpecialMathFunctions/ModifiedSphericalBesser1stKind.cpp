// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later
#include <cmath>
#include "SpecialMathFunctions/ModifiedSphericalBesser1stKind.h"

ModifiedSphericalBessel1stKind::ModifiedSphericalBessel1stKind(int degree) :
        _degree(degree) {
    _in.reserve(degree);
    _din.reserve(degree);
}

void ModifiedSphericalBessel1stKind::evaluate(double r, bool differentiate) {
    _in.clear();
    _din.clear();

    _in = ModifiedSphericalBessel1stKind::eval(_degree, r);

    if (differentiate) {
        _din.resize(_degree+1, 0.);
        if (r < RADZERO) {
            _din[1] = 1./3.;
            //_din.push_back(0.0);
            //_din.push_back(1./3.);
            //for (int n = 2; n <= _degree; ++n) {
            //    _din.push_back(0.);
            //}
        }
        else {
            _din[0] = _in[1];
            for (int n = 1; n <= _degree; ++n) {
                _din[n] = _in[n-1] - (n+1.)/r*_in[n];
            }
            //_din.push_back( _in[1] );
            //for (int n = 1; n <= _degree; ++n) {
            //    _din.push_back( _in[n-1] - (n+1.)/r*_in[n] );
            //}
        }
    }
}

std::vector<double> ModifiedSphericalBessel1stKind::eval(int degree, double r) {
    std::vector<double> il;
    if (r < RADZERO) {
        il.push_back(1.);
        il.push_back(0.);
    }
    else {
        il.push_back(sinh(r)/r);
        il.push_back(cosh(r)/r - sinh(r)/(r*r));
    }
    for (int l = 2; l <= degree; ++l) {
        if (r < RADZERO) {
            il.push_back(0.);
        }
        else {
            if (il[l-1] < SPHZERO) il.push_back(0.);
            il.push_back( il[l-2] - (2*(l-1)+1)/r*il[l-1] );
        }
    }
    return il;
}