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
    bool operator==(const ExpansionSettings& other) const;
    void checkBounds(unsigned n, unsigned l, int m) const;

    class RadialGaussianBasisSettings{
    public:
        RadialGaussianBasisSettings() = default;
        static RadialGaussianBasisSettings defaults();
        bool operator==(const RadialGaussianBasisSettings& other) const;
        void checkBounds(unsigned n) const;


        unsigned nmax;
        RadialGaussianBasisType basisType;
        double sigmaAtom, cutoffRadius;
    };

    class AngularBasisSettings{
    public:
        AngularBasisSettings() = default;
        static AngularBasisSettings defaults();
        bool operator==(const AngularBasisSettings& other) const;
        void checkBounds(unsigned l, int m) const;

        unsigned lmax;
    };

    RadialGaussianBasisSettings radial;
    AngularBasisSettings angular;
};


#endif //AMOLQCPP_EXPANSIONSETTINGS_H
