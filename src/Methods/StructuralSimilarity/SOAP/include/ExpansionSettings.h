//
// Created by Michael Heuer on 03.05.18.
//

#ifndef INPSIGHTS_EXPANSIONSETTINGS_H
#define INPSIGHTS_EXPANSIONSETTINGS_H

#include <cassert>
#include <ElementType.h>
#include <SpinType.h>
#include <NaturalConstants.h>

#include <ISettings.h>

namespace Settings{
    namespace Alchemical{
        extern std::map<std::pair<int,int>,double> pairSimilarities;
        std::string toString();
    }
};

namespace SOAPExpansion {
    enum class Mode {
        undefined = -1,
        typeAgnostic = 0, chemical = 1, alchemical = 2
    };
}

namespace Settings{
    class SOAPExpansion : public ISettings {
        inline static const std::string className = {VARNAME(SOAPExpansion)};
    public:
        Property<::SOAPExpansion::Mode> mode = {::SOAPExpansion::Mode::typeAgnostic, VARNAME(mode)};
        Property<double> zeta = {2.0, VARNAME(zeta)};
        Property<double> gamma = {0.1, VARNAME(gamma)};

        SOAPExpansion() = default;
        explicit SOAPExpansion(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };

    class Angular : public ISettings {
        inline static const std::string className = {VARNAME(Angular)};
    public:
        Property<unsigned> lmax = {5, VARNAME(lmax)};

        Angular() = default;
        explicit Angular(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };


    class Radial : public ISettings {
        inline static const std::string className = {VARNAME(Radial)};
    public:
        Property<double> radiusZero = 1e-10;

        Property<unsigned> nmax = {5, VARNAME(nmax)};
        //extern BasisType basisType;
        Property<double> sigmaAtom = {0.5*ConversionFactors::angstrom2bohr, VARNAME(sigmaAtom)};
        Property<unsigned> integrationSteps = {100, VARNAME(integrationSteps)};
        Property<double> desiredAbsoluteError = {0.0, VARNAME(desiredAbsoluteError)};
        Property<double> desiredRelativeError = {std::numeric_limits<double>::epsilon()*1e2, VARNAME(desiredRelativeError)};

        Radial() = default;
        explicit Radial(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };


    class Cutoff : public ISettings {
        inline static const std::string className = {VARNAME(Cutoff)};
    public:
        Property<double> radius = {4.0*ConversionFactors::angstrom2bohr, VARNAME(radius)};
        Property<double> width = {1.0*ConversionFactors::angstrom2bohr, VARNAME(width)};
        Property<double> centerWeight = {1.0, VARNAME(centerWeight)}; //TODO

        Cutoff() = default;
        explicit Cutoff(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::Angular)
YAML_SETTINGS_DECLARATION(Settings::Radial)
YAML_SETTINGS_DECLARATION(Settings::Cutoff)
YAML_SETTINGS_DECLARATION(Settings::SOAPExpansion)


namespace SOAPExpansion{
    extern Settings::SOAPExpansion settings;

    void checkBounds(unsigned n, unsigned l, int m);

    std::string toString(Mode mode);

    Mode fromString(std::string string);
}

namespace Angular {
    extern Settings::Angular settings;
    void checkBounds(unsigned l, int m = 0);
};

namespace Radial {
    extern Settings::Radial settings;

    enum class BasisType {
        equispaced = 0, adaptive = 1,
    };
    extern BasisType basisType;

    std::string toString(const BasisType &type);

    std::ostream &operator<<(std::ostream &os, const BasisType &type);

    void checkBounds(unsigned n);
};

namespace Cutoff{
    extern Settings::Cutoff settings;
    double innerPlateauRadius();
}



#endif //INPSIGHTS_EXPANSIONSETTINGS_H
