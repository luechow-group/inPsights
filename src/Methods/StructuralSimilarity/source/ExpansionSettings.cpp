//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>
#include <limits>
#include <SpinType.h>
#include <sstream>
#include <ElementInfo.h>

namespace ExpansionSettings {

    ExpansionSettings::Mode mode = ExpansionSettings::Mode::chemical;
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
        mode = ExpansionSettings::Mode::chemical;
        zeta = 2;
        gamma = 1.0; //TODO find sensible default value

    };

    std::string toString(const Mode &mode) {
        switch(mode) {
            case ExpansionSettings::Mode::generic : return "generic";
            case ExpansionSettings::Mode::chemical : return "chemical";
            case ExpansionSettings::Mode::alchemical : return "alchemical";
        }
    }

    std::ostream &operator<<(std::ostream &os, const Mode &mode) {
        os << toString(mode);
        return os;
    }

    std::string toString() {
        std::stringstream ss;
        ss << "General:" << std::endl
           << "--------" << std::endl
           << "Expansion mode\t\t: "<< ExpansionSettings::mode << std::endl
           << "Sharpness zeta\t\t: "<< zeta << std::endl
           << "Regularization gamma: "<< gamma << std::endl
           << std::endl
           << Radial::toString() << std::endl
           << Angular::toString() << std::endl
           << Cutoff::toString() << std::endl
           << Alchemical::toString() << std::endl;
        return ss.str();
    }

    namespace Radial {
        unsigned nmax = 14;
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

        std::string toString(const BasisType &type) {
            switch(type) {
                case BasisType::equispaced : return "equispaced";
                case BasisType::adaptive : return "adaptive";
            }
        }

        std::ostream &operator<<(std::ostream &os, const BasisType &type) {
            os << toString(type);
            return os;
        }

        std::string toString() {
            std::stringstream ss;
            ss << "Radial:" << std::endl
               << "-------" << std::endl
               << "BasisType\t\t\t: " << basisType << std::endl
               << "n_max\t\t\t\t: " << nmax << std::endl
               << "sigma_atom\t\t\t: " << sigmaAtom << " angstrom" << std::endl
               << "Integration steps\t: " << integrationSteps << std::endl
               << "Desired abs. err\t: " << desiredAbsoluteError << std::endl
               << "Desired rel. err\t: " << desiredRelativeError << std::endl;
            return ss.str();
        }
    }

    namespace Angular {
        unsigned lmax = 14;

        void defaults() {
            lmax = 14;
        };

        void checkBounds(unsigned l, int m) {
            assert(l <= lmax && "l must be less than or equal to lmax");
            assert(unsigned(abs(m)) <= lmax && "abs(m) must be smaller than lmax");
        }

        std::string toString() {
            std::stringstream ss;
            ss << "Angular:" << std::endl
               << "--------" << std::endl
               << "l_max\t\t\t\t: " << lmax  << std::endl;
            return ss.str();
        }
    }
    
    namespace Cutoff {
        double radius = 8.0;//4.0*ConversionFactors::angstrom2bohr;
        double width = 2.0;//1.0*ConversionFactors::angstrom2bohr; //TODO
        double centerWeight = 1.0; //TODO

        void defaults() {
            radius = 8.0;//4.0*ConversionFactors::angstrom2bohr;
            width = 2.0;//1.0*ConversionFactors::angstrom2bohr;
            centerWeight = 1.0;
        }

        double innerPlateauRadius() {
            return radius - width;
        }

        std::string toString() {
            std::stringstream ss;
            ss << "Cutoff:" << std::endl
               << "-------" << std::endl
               << "Radius\t\t\t\t: " << radius << " angstrom" << std::endl
               << "Width\t\t\t\t: " << width << " angstrom" << std::endl
               << "Center weight\t\t: " << centerWeight << std::endl;
            return ss.str();
        }
    }
    namespace Alchemical{
        std::map<std::pair<int,int>,double> pairSimilarities = {
                {{int(Spin::alpha),int(Spin::beta)}, 0.5} //pairs start with the smaller typeId
        };

        std::string toString() {
            std::stringstream ss;
            ss << "Alchemical Similarities:" << std::endl
               << "------------------------" << std::endl;

            std::map<std::pair<int,int>, double>::iterator it;

            for (it = pairSimilarities.begin(); it != pairSimilarities.end(); it++)
            {
                int type1 = it->first.first;
                if(type1 < 0) ss << Spins::toString(Spin(type1));
                else ss << Elements::ElementInfo::symbol(Element(type1));

                ss << " <-> ";

                int type2 = it->first.second;
                if(type2 < 0) ss << Spins::toString(Spin(type2));
                else ss << Elements::ElementInfo::symbol(Element(type2));

                ss << ": " << it->second << std::endl;
            }

            return ss.str();
        }
    }
}
