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
    SOAPExpansion::SOAPExpansion(const YAML::Node &node) {
        auto modeNode = node[className][mode.name()];
        if(!modeNode)
            spdlog::info("Property \"{0}\" was not found. Using preset value: {1}",
                    mode.name(), ::SOAPExpansion::toString(mode.get()));
        else
            mode = ::SOAPExpansion::fromString(modeNode.as<std::string>());

        auto mapNode = node[className][VARNAME(pairSimilarities)];
        if(!mapNode){
            spdlog::info("Property \"{0}\" was not found. Using preset values:");
            for(const auto& kv : pairSimilarities) {
                if (kv.first.first < 0)
                    spdlog::info("\t{0},{1} : {2}", Spins::toString(Spins::SpinType(kv.first.first)));
                else
                    spdlog::info("\t{0},{1} : {2}", Elements::ElementInfo::symbol(Elements::ElementType(kv.first.first)));

                if (kv.first.second < 0)
                    spdlog::info("\t{0},{1} : {2}", Spins::toString(Spins::SpinType(kv.first.second)));
                else
                    spdlog::info("\t{0},{1} : {2}", Elements::ElementInfo::symbol(Elements::ElementType(kv.first.second)));
            }


        } else {
            pairSimilarities.clear();
            for (const auto &kv : mapNode) {
                auto symbolPair = kv.first.as<std::vector<std::string>>();

                int type1;
                if (symbolPair[0] == Spins::toString(Spins::SpinType::alpha))
                    type1 = int(Spins::SpinType::alpha);
                else if (symbolPair[0] == Spins::toString(Spins::SpinType::beta))
                    type1 = int(Spins::SpinType::beta);
                else
                    type1 = int(Elements::ElementInfo::elementTypeFromSymbol(symbolPair[0]));

                int type2;
                if (symbolPair[1] == Spins::toString(Spins::SpinType::alpha))
                    type2 = int(Spins::SpinType::alpha);
                else if (symbolPair[1] == Spins::toString(Spins::SpinType::beta))
                    type2 = int(Spins::SpinType::beta);
                else
                    type2 = int(Elements::ElementInfo::elementTypeFromSymbol(symbolPair[1]));

                pairSimilarities.emplace(std::pair<int, int>(type1, type2), kv.second.as<double>());
            }
        }


        doubleProperty::decode(node[className], zeta);
        doubleProperty::decode(node[className], gamma);
    }

    void SOAPExpansion::appendToNode(YAML::Node &node) const {
        YAML::Node mapNode;

        for(const auto& pair : pairSimilarities) {
            std::vector<std::string> symbolPair;

            if(pair.first.first < 0)
                symbolPair.emplace_back(Spins::toString(Spins::SpinType (pair.first.first)));
            else
                symbolPair.emplace_back(Elements::ElementInfo::symbol(Elements::ElementType(pair.first.first)));

            if(pair.first.second < 0)
                symbolPair.emplace_back(Spins::toString(Spins::SpinType (pair.first.second)));
            else
                symbolPair.emplace_back(Elements::ElementInfo::symbol(Elements::ElementType(pair.first.second)));

            mapNode[symbolPair] = pair.second;
        }

        node[className][VARNAME(pairSimilarities)] = mapNode;
        node[className][mode.name()] = ::SOAPExpansion::toString(::SOAPExpansion::Mode(mode.get()));
        node[className][zeta.name()] = zeta.get();
        node[className][gamma.name()] = gamma.get();
    }


    Angular::Angular(const YAML::Node &node) {
        unsignedProperty::decode(node[className], lmax);
    }

    void Angular::appendToNode(YAML::Node &node) const {
        node[className][lmax.name()] = lmax.get();
    }


    Radial::Radial(const YAML::Node &node) {
        auto basisTypeNode = node[className][basisType.name()];
        if(!basisTypeNode)
            spdlog::info("Property \"{0}\" was not found. Using preset value: {1}",
                         basisType.name(), ::Radial::toString(basisType.get()));
        else
            basisType = ::Radial::fromString(basisTypeNode.as<std::string>());

        unsignedProperty::decode(node[className], nmax);
        doubleProperty::decode(node[className], sigmaAtom);
        unsignedProperty::decode(node[className], integrationSteps);
        doubleProperty::decode(node[className], desiredAbsoluteError);
        doubleProperty::decode(node[className], desiredRelativeError);
    }

    void Radial::appendToNode(YAML::Node &node) const {
        node[className][basisType.name()] = ::Radial::toString(::Radial::BasisType(basisType.get()));
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
YAML_SETTINGS_DEFINITION(Settings::SOAPExpansion)

namespace SOAPExpansion{
    Settings::SOAPExpansion settings = Settings::SOAPExpansion();

    std::string toString(SOAPExpansion::Mode mode) {
        switch(mode) {
            case SOAPExpansion::Mode::chemical :
                return "chemical";
            case SOAPExpansion::Mode::alchemical :
                return "alchemical";
            case SOAPExpansion::Mode::typeAgnostic :
                return "typeAgnostic";
            default:
                return "undefined";
        }
    }

    Mode fromString(const std::string& string) {
        if(string == "chemical")
            return SOAPExpansion::Mode::chemical;
        else if (string == "alchemical")
            return SOAPExpansion::Mode::alchemical;
        else if (string == "typeAgnostic")
            return SOAPExpansion::Mode::typeAgnostic;
        else
            return SOAPExpansion::Mode::undefined;
    }

    void checkBounds(unsigned n, unsigned l, int m) {
        ::Radial::checkBounds(n);
        ::Angular::checkBounds(l, m);
    }
}

namespace Angular{
    Settings::Angular settings = Settings::Angular();

    void checkBounds(unsigned l, int m) {
        assert(l <= settings.lmax.get() && "l must be less than or equal to lmax");
        assert(unsigned(abs(m)) <= settings.lmax.get() && "abs(m) must be smaller than lmax");
    }
}

namespace Radial{
    Settings::Radial settings = Settings::Radial();

    void checkBounds(unsigned n) {
        assert(n <= settings.nmax.get() && "n must be smaller than nmax");
        assert(n >= 1 && "n must greater than or equal to 1");
    }

    std::string toString(BasisType type) {
        switch(type) {
            case BasisType::equispaced :
                return "equispaced";
            case BasisType::adaptive :
                return "adaptive";
            default:
                return "undefined";
        }
    }

    BasisType fromString(const std::string& string) {
        if(string == "equispaced")
            return Radial::BasisType::equispaced;
        else if (string == "adaptive")
            return Radial::BasisType::adaptive;
        else
            return Radial::BasisType::undefined;
    }
}

namespace Cutoff{
    Settings::Cutoff settings = Settings::Cutoff();

    double innerPlateauRadius() {
        return settings.radius.get() - settings.width.get();
    }
}
