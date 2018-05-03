//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_EXPANSIONSETTINGS_H
#define AMOLQCPP_EXPANSIONSETTINGS_H

enum class RadialGaussianBasisType{
    equispaced = 0, adaptive,
};


class RadialGaussianBasisSettings{
public:
    RadialGaussianBasisSettings() = default;

    static RadialGaussianBasisSettings defaults(){
        RadialGaussianBasisSettings s{};
        s.nmax = 10;
        s.basisType = RadialGaussianBasisType::equispaced;
        s.sigmaAtom = 0.5;
        s.cutoffRadius = 4.0;

        return s;
    };

    unsigned nmax;
    RadialGaussianBasisType basisType;
    double sigmaAtom, cutoffRadius;
};

class AngularBasisSettings{
public:
    AngularBasisSettings() = default;

    static AngularBasisSettings defaults(){
        AngularBasisSettings s{};
        s.lmax = 10;
        return s;
    };

    unsigned lmax;
};

class ExpansionSettings{
public:
    ExpansionSettings() = default;

    static ExpansionSettings defaults(){
        ExpansionSettings s{};
        s.radial = RadialGaussianBasisSettings::defaults();
        s.angular = AngularBasisSettings::defaults();

        return s;
    };

    RadialGaussianBasisSettings radial;
    AngularBasisSettings angular;
};


#endif //AMOLQCPP_EXPANSIONSETTINGS_H
