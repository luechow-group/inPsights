//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_EXPANSIONSETTINGS_H
#define AMOLQCPP_EXPANSIONSETTINGS_H

enum class RadialGaussianBasisType{
    equispaced = 0, adaptive,
};

namespace ZeroLimits{
    const double radiusZero = 1e-10; //TODO put in expansion settings
}

class ExpansionSettings{
public:
    static ExpansionSettings defaults();
    bool operator==(const ExpansionSettings& other) const;
    void checkBounds(unsigned n, unsigned l, int m) const;

    class RadialGaussianBasisSettings{
    public:
        static RadialGaussianBasisSettings defaults();
        bool operator==(const RadialGaussianBasisSettings& other) const;
        void checkBounds(unsigned n) const;

        unsigned nmax;
        RadialGaussianBasisType basisType;
        double sigmaAtom, cutoffRadius;

    private:
        RadialGaussianBasisSettings() = default;
    };

    class AngularBasisSettings{
    public:
        static AngularBasisSettings defaults();
        bool operator==(const AngularBasisSettings& other) const;
        void checkBounds(unsigned l, int m) const;

        unsigned lmax;

    private:
        AngularBasisSettings() = default;
    };

    RadialGaussianBasisSettings radial;
    AngularBasisSettings angular;

private:
    ExpansionSettings() = default;
};


#endif //AMOLQCPP_EXPANSIONSETTINGS_H
