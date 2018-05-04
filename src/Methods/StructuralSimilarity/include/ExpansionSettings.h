//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_EXPANSIONSETTINGS_H
#define AMOLQCPP_EXPANSIONSETTINGS_H

enum class RadialGaussianBasisType{
    equispaced = 0, adaptive,
};


class ExpansionSettings{
public:
    ExpansionSettings() = default;

    static ExpansionSettings defaults();

    class RadialGaussianBasisSettings{
    public:
        RadialGaussianBasisSettings() = default;
        static RadialGaussianBasisSettings defaults();

        unsigned nmax;
        RadialGaussianBasisType basisType;
        double sigmaAtom, cutoffRadius;
    };

    class AngularBasisSettings{
    public:
        AngularBasisSettings() = default;
        static AngularBasisSettings defaults();

        unsigned lmax;
    };

    RadialGaussianBasisSettings radial;
    AngularBasisSettings angular;
};


#endif //AMOLQCPP_EXPANSIONSETTINGS_H
