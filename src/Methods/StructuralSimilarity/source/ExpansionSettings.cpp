//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>

namespace ExpansionSettings {
    unsigned Radial::nmax = 0;
    RadialGaussianBasisType Radial::basisType = RadialGaussianBasisType::equispaced;
    double Radial::sigmaAtom = 0;

    unsigned Angular::lmax = 0;

    double Cutoff::cutoffRadius = 0;
    double Cutoff::cutoffWidth = 0;
    double Cutoff::centerWeight = 0;

    ExpansionMode mode = ExpansionMode::Generic;


    void defaults() {
        Radial::defaults();
        Angular::defaults();
        Cutoff::defaults();
        mode = ExpansionMode::Generic;
    };

    void Radial::defaults() {
        nmax = 5;
        basisType = RadialGaussianBasisType::equispaced;
        sigmaAtom = 0.5;
    };

    void Angular::defaults() {
        lmax = 1;
    };

    void Cutoff::defaults() {
        cutoffRadius = 4.0;
        cutoffWidth = 1.0;
        centerWeight = 1.0;
    }

    double Cutoff::innerPlateauRadius() {
        return cutoffRadius - cutoffWidth;
    }

    void checkBounds(unsigned n, unsigned l, int m) {
        Radial::checkBounds(n);
        Angular::checkBounds(l,m);
    }

    void Radial::checkBounds(unsigned n) {
        assert( n <= nmax && "n must be smaller than nmax");
        assert( n >= 1 && "n must greater than or equal to 1");
    }

    void Angular::checkBounds(unsigned l, int m) {
        assert( l <= lmax && "l must be less than or equal to lmax");
        assert( unsigned(abs(m)) <= lmax && "abs(m) must be smaller than lmax");
    }
}
