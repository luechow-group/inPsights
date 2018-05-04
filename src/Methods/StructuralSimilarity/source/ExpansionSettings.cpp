//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>

ExpansionSettings ExpansionSettings::defaults() {
    ExpansionSettings s{};
    s.radial = ExpansionSettings::RadialGaussianBasisSettings::defaults();
    s.angular = ExpansionSettings::AngularBasisSettings::defaults();

    return s;
};

ExpansionSettings::RadialGaussianBasisSettings ExpansionSettings::RadialGaussianBasisSettings::defaults() {
    ExpansionSettings::RadialGaussianBasisSettings s{};
    s.nmax = 10;
    s.basisType = RadialGaussianBasisType::equispaced;
    s.sigmaAtom = 0.5;
    s.cutoffRadius = 4.0;

    return s;
};

ExpansionSettings::AngularBasisSettings ExpansionSettings::AngularBasisSettings::defaults() {
    ExpansionSettings::AngularBasisSettings s{};
    s.lmax = 10;

    return s;
};

bool ExpansionSettings::operator==(const ExpansionSettings &other) const {
    return (this->radial == other.radial) && (this->angular == other.angular);
}

bool ExpansionSettings::RadialGaussianBasisSettings::operator==(
        const ExpansionSettings::RadialGaussianBasisSettings& other) const {
    return (this->nmax == other.nmax) &&
           (this->basisType == other.basisType) &&
           (this->sigmaAtom == other.sigmaAtom) &&
           (this->cutoffRadius == other.cutoffRadius);
}

bool ExpansionSettings::AngularBasisSettings::operator==(
        const ExpansionSettings::AngularBasisSettings& other) const {
    return this->lmax == other.lmax;
}

void ExpansionSettings::checkBounds(unsigned n, unsigned l, int m) const {
    radial.checkBounds(n);
    angular.checkBounds(l,m);
}

void ExpansionSettings::RadialGaussianBasisSettings::checkBounds(unsigned n) const {
    assert( n <= nmax && "n must be smaller than nmax");
    assert( n >= 1 && "n must greater than or equal to 1");
}

void ExpansionSettings::AngularBasisSettings::checkBounds(unsigned l, int m) const {
    assert( l <= lmax && "l must be smaller than lmax");
    assert( unsigned(abs(m)) <= lmax && "abs(m) must be smaller than lmax");
}
