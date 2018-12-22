//
// Created by Michael Heuer on 04.05.18.
//

#include "ExpansionSettings.h"
#include <cassert>
#include <limits>
#include <SpinType.h>
#include <sstream>
#include <ElementInfo.h>
#include <ExpansionSettings.h>


namespace ExpansionSettings {
    void checkBounds(unsigned n, unsigned l, int m) {
        Radial::checkBounds(n);
        Angular::checkBounds(l, m);
    }

    std::string toString(const Mode &mode) {
        switch(mode) {
            case ExpansionSettings::Mode::typeAgnostic :
                return "typeAgnostic";
            case ExpansionSettings::Mode::chemical :
                return "chemical";
            case ExpansionSettings::Mode::alchemical :
                return "alchemical";
            default:
                return "undefined";
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
        void checkBounds(unsigned n) {
            assert(n <= nmax && "n must be smaller than nmax");
            assert(n >= 1 && "n must greater than or equal to 1");
        }

        std::string toString(const BasisType &type) {
            switch(type) {
                case BasisType::equispaced :
                    return "equispaced";
                case BasisType::adaptive :
                    return "adaptive";
                default:
                    return "undefined";
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
               << "sigma_atom\t\t\t: " << sigmaAtom << " bohr" << std::endl
               << "Integration steps\t: " << integrationSteps << std::endl
               << "Desired abs. err\t: " << desiredAbsoluteError << std::endl
               << "Desired rel. err\t: " << desiredRelativeError << std::endl;
            return ss.str();
        }
    }

    namespace Angular {

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
        double innerPlateauRadius() {
            return radius - width;
        }

        std::string toString() {
            std::stringstream ss;
            ss << "Cutoff:" << std::endl
               << "-------" << std::endl
               << "Radius\t\t\t\t: " << radius << " bohr" << std::endl
               << "Width\t\t\t\t: " << width << " bohr" << std::endl
               << "Center weight\t\t: " << centerWeight << std::endl;
            return ss.str();
        }
    }
    namespace Alchemical{
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
