//
// Created by Michael Heuer on 03.05.18.
//

#ifndef INPSIGHTS_EXPANSIONSETTINGS_H
#define INPSIGHTS_EXPANSIONSETTINGS_H

#include <cassert>
#include <ElementType.h>
#include <SpinType.h>
#include <NaturalConstants.h>

namespace ExpansionSettings{
    namespace Radial{
        enum class BasisType{
            equispaced = 0, adaptive,
        };

        std::string toString(const BasisType& type);
        std::ostream& operator<<(std::ostream& os, const BasisType& type);

        const double radiusZero = 1e-10;

        inline unsigned nmax = {5};
        inline BasisType basisType = {BasisType::equispaced};
        inline double sigmaAtom = {0.5*ConversionFactors::angstrom2bohr};

        inline unsigned integrationSteps = {100};
        inline double desiredAbsoluteError = {0.0}, desiredRelativeError = {std::numeric_limits<double>::epsilon()*1e2};

        void checkBounds(unsigned n);

        std::string toString();
    };

    namespace Angular {
        inline unsigned lmax = {5};

        void checkBounds(unsigned l, int m = 0);

        std::string toString();
    };

    namespace Cutoff {
        inline double radius = {4.0*ConversionFactors::angstrom2bohr};
        inline double width = {1.0*ConversionFactors::angstrom2bohr};
        inline double centerWeight = {1.0}; //TODO

        double innerPlateauRadius();

        std::string toString();
    }

    namespace Alchemical{
        inline std::map<std::pair<int,int>,double> pairSimilarities = {
                {{int(Spin::alpha),int(Spin::beta)}, 0.5}
        };
        std::string toString();
    }

    enum class Mode {
        typeAgnostic = 0, chemical, alchemical
    };

    inline Mode mode = ExpansionSettings::Mode::chemical;
    inline double zeta = 2; // LocalSimilarity exponent
    inline double gamma = 0.1; // StructuralSimilarity regularization parameter

    void checkBounds(unsigned n, unsigned l, int m);

    std::string toString(const Mode& mode);

    std::ostream& operator<<(std::ostream& os, const Mode& mode);

    std::string toString();
};

#endif //INPSIGHTS_EXPANSIONSETTINGS_H
