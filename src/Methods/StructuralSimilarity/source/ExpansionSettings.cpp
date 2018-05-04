//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"


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

