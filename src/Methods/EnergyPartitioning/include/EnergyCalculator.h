//
// Created by Michael Heuer on 02.10.18.
//

#ifndef INPSIGHTS_ENERGYCALCULATOR_H
#define INPSIGHTS_ENERGYCALCULATOR_H

#include <utility>
#include <Statistics.h>
#include <Reference.h>
#include <Sample.h>
#include <Logger.h>
#include <SimilarReferences.h>

class EnergyCalculator{
public:

    struct TotalEnergies{
        TotalEnergies()
        : Te(0),Vee(0),Ven(0),Vnn(0){};

        double totalEnergy() const {
            return Te + Vee + Ven + Vnn;
        };

        double Te, Vee, Ven, Vnn;
    };

    EnergyCalculator(const std::vector<Sample>& samples, AtomsVector atoms);

    TotalEnergies calculateTotalEnergies();

    unsigned long addReference(const Reference &reference);

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima);

    // selects nWanted structures and prints the statistic data
    void printCluster(std::vector<ElectronsVector>& structures);

    YAML::Node getYamlNode();

    std::string getYamlDocumentString();

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    SingleParticlesStatistics valueStats_,TeStats_;
    IntraParticlesStatistics SeeStats_, VeeStats_, VnnStats_;
    InterParticlesStatistics VenStats_;

    Eigen::MatrixXd Vnn_;

    YAML::Emitter yamlDocument_;
    std::shared_ptr<spdlog::logger> console;
};

#endif //INPSIGHTS_ENERGYCALCULATOR_H
