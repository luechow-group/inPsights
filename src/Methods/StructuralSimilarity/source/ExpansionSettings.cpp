//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>

// undefined settings

unsigned ExpansionSettings::Radial::nmax = 0;
RadialGaussianBasisType ExpansionSettings::Radial::basisType = RadialGaussianBasisType(0);
double ExpansionSettings::Radial::sigmaAtom = 0;
double ExpansionSettings::Radial::cutoffRadius = 0;
unsigned ExpansionSettings::Angular::lmax = 0;

double ExpansionSettings::Cutoff::cutoffRadius = 0;
double ExpansionSettings::Cutoff::cutoffWidth = 0;
double ExpansionSettings::Cutoff::centerWeight = 0;

ExpansionMode ExpansionSettings::mode = ExpansionMode::Generic;


ExpansionSettings ExpansionSettings::defaults() {
    ExpansionSettings s{};

    s.radial = ExpansionSettings::Radial::defaults();
    s.angular = ExpansionSettings::Angular::defaults();
    s.cutoff = ExpansionSettings::Cutoff::defaults();
    ExpansionSettings::mode = ExpansionMode::Generic;

    return s;
};

ExpansionSettings::Radial
ExpansionSettings::Radial::defaults() {
    ExpansionSettings::Radial r{};

    ExpansionSettings::Radial::nmax = 2;
    ExpansionSettings::Radial::basisType = RadialGaussianBasisType::equispaced;
    ExpansionSettings::Radial::sigmaAtom = 0.5;
    ExpansionSettings::Radial::cutoffRadius = 4.0;

    return r;
};

ExpansionSettings::Angular
ExpansionSettings::Angular::defaults() {
    ExpansionSettings::Angular a{};
    ExpansionSettings::Angular::lmax = 1;

    return a;
};

ExpansionSettings::Cutoff
ExpansionSettings::Cutoff::defaults() {
    ExpansionSettings::Cutoff c{};

    ExpansionSettings::Cutoff::cutoffRadius = 4.0;
    ExpansionSettings::Cutoff::cutoffWidth = 1.0;
    ExpansionSettings::Cutoff::centerWeight = 1.0;

    return c;
}

void ExpansionSettings::checkBounds(unsigned n, unsigned l, int m) {
    ExpansionSettings::Radial::checkBounds(n);
    ExpansionSettings::Angular::checkBounds(l,m);
}

void ExpansionSettings::Radial::checkBounds(unsigned n) {
    assert( n <= nmax && "n must be smaller than nmax");
    assert( n >= 1 && "n must greater than or equal to 1");
}

void ExpansionSettings::Angular::checkBounds(unsigned l, int m) {
    assert( l <= lmax && "l must be less than or equal to lmax");
    assert( unsigned(abs(m)) <= lmax && "abs(m) must be smaller than lmax");
}
