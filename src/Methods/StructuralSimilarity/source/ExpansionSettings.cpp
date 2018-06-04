//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>
#include <limits>
#include <SpinType.h>

namespace ExpansionSettings {

    ExpansionSettings::Mode mode = ExpansionSettings::Mode::TypeSpecific;
    double zeta = 2;
    double gamma = 1.0;

    void checkBounds(unsigned n, unsigned l, int m) {
        Radial::checkBounds(n);
        Angular::checkBounds(l, m);
    }

    void defaults() {
        Radial::defaults();
        Angular::defaults();
        Cutoff::defaults();
        mode = ExpansionSettings::Mode::TypeSpecific;
        zeta = 2;
        gamma = 1.0; //TODO find sensible default value

    };

    namespace Radial {
        unsigned nmax = 5;
        BasisType basisType = BasisType::equispaced;
        double sigmaAtom = 0.5;

        unsigned integrationSteps = 100;
        double desiredAbsoluteError = 0.0;
        double desiredRelativeError = 1e-6; //TODO

        void defaults() {
            nmax = 5;
            basisType = BasisType::equispaced;
            sigmaAtom = 0.5;
            integrationSteps = 100;
            desiredAbsoluteError = 0.0;
            desiredRelativeError = std::numeric_limits<double>::epsilon()*1e2;
        };

        void checkBounds(unsigned n) {
            assert(n <= nmax && "n must be smaller than nmax");
            assert(n >= 1 && "n must greater than or equal to 1");
        }

    }

    namespace Angular {
        unsigned lmax = 3;

        void defaults() {
            lmax = 3;
        };

        void checkBounds(unsigned l, int m) {
            assert(l <= lmax && "l must be less than or equal to lmax");
            assert(unsigned(abs(m)) <= lmax && "abs(m) must be smaller than lmax");
        }
    }
    
    namespace Cutoff {
        double cutoffRadius = 4.0;
        double cutoffWidth = 1.0; //TODO
        double centerWeight = 1.0; //TODO

        void defaults() {
            cutoffRadius = 4.0;
            cutoffWidth = 1.0;
            centerWeight = 1.0;
        }

        double innerPlateauRadius() {
            return cutoffRadius - cutoffWidth;
        }

        //namespace Alchemical{
        //    std::map<std::pair<int,int>,double> pairSimilarities = {
        //            {{int(Spin::alpha),int(Spin::beta)}, 0.5} //pairs start with the smaller typeId
        //    };
        //}
    }
}
