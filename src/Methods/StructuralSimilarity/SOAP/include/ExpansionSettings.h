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

        extern unsigned nmax;
        extern BasisType basisType;
        extern double sigmaAtom;

        extern unsigned integrationSteps;
        extern double desiredAbsoluteError;

        void checkBounds(unsigned n);

        std::string toString();
    };

    namespace Angular {
        extern unsigned lmax;

        void checkBounds(unsigned l, int m = 0);

        std::string toString();
    };

    namespace Cutoff {
        extern double radius;
        extern double width;
        extern double centerWeight; //TODO

        double innerPlateauRadius();

        std::string toString();
    }

    namespace Alchemical{
        extern std::map<std::pair<int,int>,double> pairSimilarities;
        std::string toString();
    }

    enum class Mode {
        typeAgnostic = 0, chemical, alchemical
    };

    extern Mode mode;
    extern double zeta;
    extern double gamma;

    void checkBounds(unsigned n, unsigned l, int m);

    std::string toString(const Mode& mode);

    std::ostream& operator<<(std::ostream& os, const Mode& mode);

    std::string toString();
};

#endif //INPSIGHTS_EXPANSIONSETTINGS_H
