//
// Created by Michael Heuer on 02.10.18.
//

#ifndef INPSIGHTS_ENERGYCALCULATOR_H
#define INPSIGHTS_ENERGYCALCULATOR_H

#include "Reference.h"
#include "SimilarReferences.h"
#include "Sample.h"
#include <Statistics.h>

class EnergyCalculator {
public:

    EnergyCalculator(YAML::Emitter& yamlDocument, const std::vector<Sample> &samples, AtomsVector atoms);

    unsigned long addReference(const Reference &reference);

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>> &clusteredGloballySimilarMaxima);

    // selects nWanted structures and prints the statistic data
    void printCluster(std::vector<ElectronsVector> &structures);

    YAML::Node getYamlNode();

    std::string getYamlDocumentString();

private:
    YAML::Emitter& yamlDocument_;
    const std::vector<Sample> &samples_;
    AtomsVector atoms_;

    SingleParticlesStatistics valueStats_, TeStats_, EeStats_, EnStats_;
    IntraParticlesStatistics SeeStats_, VeeStats_, VnnStats_;
    InterParticlesStatistics VenStats_;

    Eigen::MatrixXd Vnn_;
};

#endif //INPSIGHTS_ENERGYCALCULATOR_H
