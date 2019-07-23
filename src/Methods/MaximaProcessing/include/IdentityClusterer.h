//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_IDENTITYCLUSTERER_H
#define INPSIGHTS_IDENTITYCLUSTERER_H

#include "Sample.h"
#include "IClusterer.h"
#include <ISettings.h>

namespace Settings {
    class IdentityClusterer : public ISettings {
    public:
        Property<double> identityRadius = {0.01, VARNAME(identityRadius)};
        Property<double> identityValueIncrement = {1e-7, VARNAME(identityValueIncrement)};

        IdentityClusterer();
        explicit IdentityClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::IdentityClusterer)

class IdentityClusterer : public IClusterer{
public:
    static Settings::IdentityClusterer settings;

    IdentityClusterer(std::vector<Sample> &samples);
    void cluster(Group& group) override;

private:
    void subLoop(Group& group,
            Group::iterator &beginIt,
            Group::iterator &it,
            Group::iterator &endIt,
            double distThresh,
            double valueIncrement);

    // TODO This method should be located inside of a reference container class
    void addReference(Group& group,
            const Group::iterator &beginIt,
            Group::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const;

    std::vector<Sample> &samples_;
};


#endif //INPSIGHTS_IDENTITYCLUSTERER_H
