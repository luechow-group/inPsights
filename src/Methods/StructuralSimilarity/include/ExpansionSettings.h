//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_EXPANSIONSETTINGS_H
#define AMOLQCPP_EXPANSIONSETTINGS_H

enum class RadialGaussianBasisType{
    equispaced = 0, adaptive,
};

enum class ExpansionMode {
    Generic = 0, TypeSpecific
};


namespace ZeroLimits{
    const double radiusZero = 1e-10; //TODO put in expansion settings
}

// Monostate pattern
class ExpansionSettings{
public:
    ExpansionSettings() = default;
    static ExpansionSettings defaults();
    static void checkBounds(unsigned n, unsigned l, int m);

    class Radial{
    public:
        Radial() = default;
        static Radial defaults();
        static void checkBounds(unsigned n);

        static unsigned nmax;
        static RadialGaussianBasisType basisType;
        static double sigmaAtom, cutoffRadius;
    };

    class Angular{
    public:
        Angular() = default;
        static Angular defaults();
        static void checkBounds(unsigned l, int m = 0);

        static unsigned lmax;
    };

    class Cutoff{
    public:
        Cutoff() = default;
        static Cutoff defaults();

        static double cutoffRadius, cutoffWidth, centerWeight;
    };

    static ExpansionMode mode;

private:
    Radial radial;
    Angular angular;
    Cutoff cutoff;


};


#endif //AMOLQCPP_EXPANSIONSETTINGS_H
