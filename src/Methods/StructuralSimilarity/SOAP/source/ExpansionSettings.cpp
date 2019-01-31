//
// Created by Michael Heuer on 04.05.18.
//

#include <ExpansionSettings.h>
#include <cassert>
#include <limits>
#include <SpinType.h>
#include <sstream>
#include <ElementInfo.h>


namespace Settings {
    Angular::Angular(const YAML::Node &node) {
        unsignedProperty::decode(node[className], lmax);
    }

    void Angular::appendToNode(YAML::Node &node) const {
        node[className][lmax.name()] = lmax.get();
    }


    Radial::Radial(const YAML::Node &node) {
        unsignedProperty::decode(node[className], nmax);
        doubleProperty::decode(node[className], sigmaAtom);
        unsignedProperty::decode(node[className], integrationSteps);
        doubleProperty::decode(node[className], desiredAbsoluteError);
        doubleProperty::decode(node[className], desiredRelativeError);
    }

    void Radial::appendToNode(YAML::Node &node) const {
        node[className][nmax.name()] = nmax.get();
        node[className][sigmaAtom.name()] = sigmaAtom.get();
        node[className][integrationSteps.name()] = integrationSteps.get();
        node[className][desiredAbsoluteError.name()] = desiredAbsoluteError.get();
        node[className][desiredRelativeError.name()] = desiredRelativeError.get();
    }


    Cutoff::Cutoff(const YAML::Node &node) {
        doubleProperty::decode(node[className], radius);
        doubleProperty::decode(node[className], width);
        doubleProperty::decode(node[className], centerWeight);
    }

    void Cutoff::appendToNode(YAML::Node &node) const {
        node[className][radius.name()] = radius.get();
        node[className][width.name()] = width.get();
        node[className][centerWeight.name()] = centerWeight.get();
    }
}
YAML_SETTINGS_DEFINITION(Settings::Angular)
YAML_SETTINGS_DEFINITION(Settings::Radial)
YAML_SETTINGS_DEFINITION(Settings::Cutoff)


namespace Angular{
    Settings::Angular settings = Settings::Angular();

    void checkBounds(unsigned l, int m) {
        assert(l <= settings.lmax.get() && "l must be less than or equal to lmax");
        assert(unsigned(abs(m)) <= settings.lmax.get() && "abs(m) must be smaller than lmax");
    }
}

namespace Radial{
    Settings::Radial settings = Settings::Radial();
    BasisType basisType = BasisType::equispaced;

    void checkBounds(unsigned n) {
        assert(n <= settings.nmax.get() && "n must be smaller than nmax");
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
}

namespace Cutoff{
    Settings::Cutoff settings = Settings::Cutoff();

    double innerPlateauRadius() {
        return settings.radius.get() - settings.width.get();
    }
}


namespace Settings {
    Mode mode = Settings::Mode::chemical;
    double zeta = 2; // LocalSimilarity exponent
    double gamma = 0.1; // StructuralSimilarity regularization parameter

    void checkBounds(unsigned n, unsigned l, int m) {
        ::Radial::checkBounds(n);
        ::Angular::checkBounds(l, m);
    }

    std::string toString(const Mode &mode) {
        switch(mode) {
            case Settings::Mode::typeAgnostic :
                return "typeAgnostic";
            case Settings::Mode::chemical :
                return "chemical";
            case Settings::Mode::alchemical :
                return "alchemical";
            default:
                return "undefined";
        }
    }

    std::ostream &operator<<(std::ostream &os, const Mode &mode) {
        os << toString(mode);
        return os;
    }

    namespace Alchemical{
        std::map<std::pair<int,int>,double> pairSimilarities = {
                {{int(Spin::alpha),int(Spin::beta)}, 0.5}
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
