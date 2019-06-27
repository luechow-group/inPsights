//
// Created by leonard on 13.05.19.
//

#ifndef INPSIGHTS_LOCALBONDSIMILARITYCLUSTERER_H
#define INPSIGHTS_LOCALBONDSIMILARITYCLUSTERER_H

#include <IClusterer.h>
#include <ISettings.h>
#include <ParticlesVector.h>
#include <Sample.h>
#include <functional>

namespace Settings {
    class LocalBondSimilarityClusterer : public ISettings {
    public:
        Property<double> similarityRadius = {0.1, VARNAME(similarityRadius)};
        Property<long> maximalCount = {2, VARNAME(maximalCount)};
        Property<std::string> distanceMode = {"minimum", VARNAME(distanceMode)};

        LocalBondSimilarityClusterer();
        explicit LocalBondSimilarityClusterer(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::LocalBondSimilarityClusterer)

class LocalBondSimilarityClusterer : public IClusterer {
public:
    static Settings::LocalBondSimilarityClusterer settings;

    LocalBondSimilarityClusterer(std::vector<Sample> &samples, AtomsVector &nuclei,
            std::vector<Eigen::Vector3d> &positions);
    // in contrast to the other clusterers, this clusterer does no permutations
    void cluster(Group& group) override;
    std::list<long> getRelevantIndices(const ElectronsVector &electronsVector);

private:
    std::vector<Sample> &samples_;
    AtomsVector nuclei_;
    std::vector<Eigen::Vector3d> positions_;
    std::function<double(const Eigen::Vector3d &, const std::vector<Eigen::Vector3d> &)> distanceFunction_;
};

#endif //INPSIGHTS_LOCALBONDSIMILARITYCLUSTERER_H
