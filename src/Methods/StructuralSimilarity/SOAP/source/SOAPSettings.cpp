//
// Created by Michael Heuer on 04.05.18.
//

#include <SOAPSettings.h>
#include <cassert>
#include <limits>
#include <SpinType.h>
#include <sstream>
#include <ElementInfo.h>

namespace Settings {
    namespace SOAP {
        General::General()
        : ISettings(VARNAME(General)) {}

        General::General(const YAML::Node &node)
        : General() {
            auto modeNode = node[className][mode.name()];
            if (!modeNode)
                spdlog::info("Property \"{0}\" was not found. Using preset value: {1}",
                             mode.name(), ::SOAP::General::toString(mode()));
            else
                mode = ::SOAP::General::fromString(modeNode.as<std::string>());

            auto mapNode = node[className][VARNAME(pairSimilarities)];
            if (!mapNode) {
                spdlog::info("Property \"{0}\" was not found. Using preset values:");
                for (const auto &kv : pairSimilarities) {
                    if (kv.first.first < 0)
                        spdlog::info("\t{0},{1} : {2}", Spins::toString(Spins::SpinType(kv.first.first)));
                    else
                        spdlog::info("\t{0},{1} : {2}",
                                     Elements::ElementInfo::symbol(Elements::ElementType(kv.first.first)));

                    if (kv.first.second < 0)
                        spdlog::info("\t{0},{1} : {2}", Spins::toString(Spins::SpinType(kv.first.second)));
                    else
                        spdlog::info("\t{0},{1} : {2}",
                                     Elements::ElementInfo::symbol(Elements::ElementType(kv.first.second)));
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
            doubleProperty::decode(node[className], sinkhornEpsilon);
        }

        void General::appendToNode(YAML::Node &node) const {
            YAML::Node mapNode;

            for (const auto &pair : pairSimilarities) {
                std::vector<std::string> symbolPair;

                if (pair.first.first < 0)
                    symbolPair.emplace_back(Spins::toString(Spins::SpinType(pair.first.first)));
                else
                    symbolPair.emplace_back(Elements::ElementInfo::symbol(Elements::ElementType(pair.first.first)));

                if (pair.first.second < 0)
                    symbolPair.emplace_back(Spins::toString(Spins::SpinType(pair.first.second)));
                else
                    symbolPair.emplace_back(Elements::ElementInfo::symbol(Elements::ElementType(pair.first.second)));

                mapNode[symbolPair] = pair.second;
            }

            node[className][VARNAME(pairSimilarities)] = mapNode;
            node[className][mode.name()] = ::SOAP::General::toString(::SOAP::General::Mode(mode()));
            node[className][zeta.name()] = zeta();
            node[className][gamma.name()] = gamma();
            node[className][sinkhornEpsilon.name()] = sinkhornEpsilon();
        }


        Angular::Angular(): ISettings(VARNAME(Angular)) {}

        Angular::Angular(const YAML::Node &node)
        : Angular() {
            unsignedProperty::decode(node[className], lmax);
        }

        void Angular::appendToNode(YAML::Node &node) const {
            node[className][lmax.name()] = lmax();
        }

        Radial::Radial()
        : ISettings(VARNAME(Radial)) {}

        Radial::Radial(const YAML::Node &node)
        : Radial() {
            auto basisTypeNode = node[className][basisType.name()];
            if (!basisTypeNode)
                spdlog::info("Property \"{0}\" was not found. Using preset value: {1}",
                             basisType.name(), ::SOAP::Radial::toString(basisType()));
            else
                basisType = ::SOAP::Radial::fromString(basisTypeNode.as<std::string>());

            unsignedProperty::decode(node[className], nmax);
            doubleProperty::decode(node[className], sigmaAtom);
            unsignedProperty::decode(node[className], integrationSteps);
            doubleProperty::decode(node[className], desiredAbsoluteError);
            doubleProperty::decode(node[className], desiredRelativeError);
        }

        void Radial::appendToNode(YAML::Node &node) const {
            node[className][basisType.name()] = ::SOAP::Radial::toString(::SOAP::Radial::BasisType(basisType()));
            node[className][nmax.name()] = nmax();
            node[className][sigmaAtom.name()] = sigmaAtom();
            node[className][integrationSteps.name()] = integrationSteps();
            node[className][desiredAbsoluteError.name()] = desiredAbsoluteError();
            node[className][desiredRelativeError.name()] = desiredRelativeError();
        }


        Cutoff::Cutoff() : ISettings(VARNAME(Cutoff)) {}

        Cutoff::Cutoff(const YAML::Node &node)
        : Cutoff() {
            doubleProperty::decode(node[className], radius);
            doubleProperty::decode(node[className], width);
            doubleProperty::decode(node[className], centerWeight);
        }

        void Cutoff::appendToNode(YAML::Node &node) const {
            node[className][radius.name()] = radius();
            node[className][width.name()] = width();
            node[className][centerWeight.name()] = centerWeight();
        }
    }
}
YAML_SETTINGS_DEFINITION(Settings::SOAP::Angular)
YAML_SETTINGS_DEFINITION(Settings::SOAP::Radial)
YAML_SETTINGS_DEFINITION(Settings::SOAP::Cutoff)
YAML_SETTINGS_DEFINITION(Settings::SOAP::General)

namespace SOAP{
    namespace General {

        Settings::SOAP::General settings = Settings::SOAP::General();

        std::string toString(Mode mode) {
            switch (mode) {
                case Mode::chemical :
                    return "chemical";
                case Mode::alchemical :
                    return "alchemical";
                case Mode::typeAgnostic :
                    return "typeAgnostic";
                default:
                    return "undefined";
            }
        }

        Mode fromString(const std::string &string) {
            if (string == "chemical")
                return Mode::chemical;
            else if (string == "alchemical")
                return Mode::alchemical;
            else if (string == "typeAgnostic")
                return Mode::typeAgnostic;
            else
                return Mode::undefined;
        }

        void checkBounds(unsigned n, unsigned l, int m) {
            ::SOAP::Radial::checkBounds(n);
            ::SOAP::Angular::checkBounds(l, m);
        }
    }

    namespace Angular{
        Settings::SOAP::Angular settings = Settings::SOAP::Angular();

        void checkBounds(unsigned l, int m) {
            assert(l <= settings.lmax() && "l must be less than or equal to lmax");
            assert(unsigned(abs(m)) <= settings.lmax() && "abs(m) must be smaller than lmax");
        }
    }

    namespace Radial{
        Settings::SOAP::Radial settings = Settings::SOAP::Radial();

        void checkBounds(unsigned n) {
            assert(n <= settings.nmax() && "n must be smaller than nmax");
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
        Settings::SOAP::Cutoff settings = Settings::SOAP::Cutoff();

        double innerPlateauRadius() {
            return settings.radius() - settings.width();
        }
    }
}