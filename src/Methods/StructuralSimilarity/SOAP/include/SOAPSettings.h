//
// Created by Michael Heuer on 03.05.18.
//

#ifndef INPSIGHTS_SOAPSETTINGS_H
#define INPSIGHTS_SOAPSETTINGS_H

#include <cassert>
#include <ElementType.h>
#include <SpinType.h>
#include <NaturalConstants.h>

#include <ISettings.h>

namespace SOAP {
    namespace General {
        enum class Mode {
            undefined = -1,
            typeAgnostic = 0, chemical = 1, alchemical = 2
        };
    }

    namespace Radial{
        enum class BasisType {
            undefined = -1,
            equispaced = 0, adaptive = 1,
        };
    }
}

namespace Settings{
    namespace SOAP {
        class General : public ISettings {
        public:
            Property<::SOAP::General::Mode> mode = {::SOAP::General::Mode::typeAgnostic, VARNAME(mode)};
            Property<double> zeta = {2.0, VARNAME(zeta)};
            Property<double> gamma = {0.1, VARNAME(gamma)};
            Property<double> sinkhornEpsilon = {1e-8, VARNAME(sinkhornEpsilon)};

            std::map<std::pair<int, int>, double> pairSimilarities = {
                    {{int(Spin::alpha), int(Spin::beta)}, 0.5}
            };

            General();

            explicit General(const YAML::Node &node);

            void appendToNode(YAML::Node &node) const override;
        };

        class Angular : public ISettings {
        public:
            Property<unsigned> lmax = {5, VARNAME(lmax)};

            Angular();

            explicit Angular(const YAML::Node &node);

            void appendToNode(YAML::Node &node) const override;
        };


        class Radial : public ISettings {
        public:
            Property<double> radiusZero = 1e-10;

            Property<::SOAP::Radial::BasisType> basisType =
                    {::SOAP::Radial::BasisType::equispaced, VARNAME(mode)};

            Property<unsigned> nmax = {5, VARNAME(nmax)};
            Property<double> sigmaAtom = {1.0, VARNAME(sigmaAtom)};
            Property<double> sigmaZeroThreshold = {1e-10, VARNAME(sigmaZeroThreshold)};
            Property<unsigned> integrationSteps = {100, VARNAME(integrationSteps)};
            Property<double> desiredAbsoluteError = {0.0, VARNAME(desiredAbsoluteError)};
            Property<double> desiredRelativeError = {1e-6, VARNAME(desiredRelativeError)};

            Radial();

            explicit Radial(const YAML::Node &node);

            void appendToNode(YAML::Node &node) const override;
        };


        class Cutoff : public ISettings {
        public:
            Property<double> radius = {4.0 * ConversionFactors::angstrom2bohr, VARNAME(radius)};
            Property<double> width = {1.0 * ConversionFactors::angstrom2bohr, VARNAME(width)};
            Property<double> centerWeight = {1.0, VARNAME(centerWeight)}; //TODO

            Cutoff();

            explicit Cutoff(const YAML::Node &node);

            void appendToNode(YAML::Node &node) const override;
        };
    }
}
YAML_SETTINGS_DECLARATION(Settings::SOAP::General)
YAML_SETTINGS_DECLARATION(Settings::SOAP::Angular)
YAML_SETTINGS_DECLARATION(Settings::SOAP::Radial)
YAML_SETTINGS_DECLARATION(Settings::SOAP::Cutoff)


namespace SOAP {
    namespace General {
        extern Settings::SOAP::General settings;

        void checkBounds(unsigned n, unsigned l, int m);

        std::string toString(Mode mode);

        Mode fromString(const std::string &string);
    }

    namespace Angular {
        extern Settings::SOAP::Angular settings;

        void checkBounds(unsigned l, int m = 0);
    };

    namespace Radial {
        extern Settings::SOAP::Radial settings;

        std::string toString(BasisType type);

        BasisType fromString(const std::string &string);

        std::ostream &operator<<(std::ostream &os, const BasisType &type);

        void checkBounds(unsigned n);
    };

    namespace Cutoff {
        extern Settings::SOAP::Cutoff settings;

        double innerPlateauRadius();
    }
}



#endif //INPSIGHTS_SOAPSETTINGS_H
