//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_GLOBALIDENTITYSORTER_H
#define INPSIGHTS_GLOBALIDENTITYSORTER_H

#include "Reference.h"
#include "Sample.h"
#include <Logger.h>
#include <ISettings.h>
#include <Property.h>
#include <utility>
#include <vector>

namespace Settings {
    class GlobalIdentitySorter : public ISettings {
        inline static const std::string className = {VARNAME(GlobalIdentitySorter)};
    public:
        Property<double> identityRadius = {0.01, VARNAME(identityRadius)};
        Property<double> valueIncrement = {1e-7, VARNAME(valueIncrement)};

        GlobalIdentitySorter();
        explicit GlobalIdentitySorter(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::GlobalIdentitySorter)

class GlobalIdentitySorter {
public:
    inline static Settings::GlobalIdentitySorter settings;

    GlobalIdentitySorter(
            std::vector<Reference> &references,
            std::vector<Sample> &samples);
    bool sort();

private:
    void subLoop(
            std::vector<Reference>::iterator &beginIt,
            std::vector<Reference>::iterator &it,
            std::vector<Reference>::iterator &endIt,
            double distThresh,
            double valueIncrement);

    // TODO This method should be located inside of a reference container class
    void addReference(
            const std::vector<Reference>::iterator &beginIt,
            std::vector<Reference>::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const;

    std::vector<Reference> &references_;
    std::vector<Sample> &samples_;
};


#endif //INPSIGHTS_GLOBALIDENTITYSORTER_H
