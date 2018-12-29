//
// Created by heuer on 12.12.18.
//
#include <EnergyCalculator.h>

#include <ClusterData.h>
#include <SpinCorrelation.h>
#include <CoulombPotential.h>
#include <OneParticleEnergies.h>
#include <Logger.h>

EnergyCalculator::EnergyCalculator(YAML::Emitter& yamlDocument, const std::vector<Sample>& samples, AtomsVector atoms)
        :
        yamlDocument_(yamlDocument),
        samples_(samples),
        atoms_(std::move(atoms)),
        Vnn_(CoulombPotential::energies<Element>(atoms_))
{

    VnnStats_.add(Vnn_);
    EnStats_.add(OneParticleEnergies::oneAtomEnergies(Vnn_));
}

unsigned long EnergyCalculator::addReference(const Reference &reference) {
    auto count = unsigned(reference.count());

    Eigen::VectorXd value(1); value[0] = reference.value();
    Eigen::MatrixXd spinCorrelations_ = SpinCorrelation::spinCorrelations(reference.maximum().typesVector()).cast<double>();
    // Maximum related statistics
    valueStats_.add(value, count);
    SeeStats_.add(spinCorrelations_,count);

    // Sample related statistics
    for (auto & id : reference.sampleIds()) {
        Eigen::VectorXd Te = samples_[id].kineticEnergies_;
        Eigen::MatrixXd Vee = CoulombPotential::energies(samples_[id].sample_);
        Eigen::MatrixXd Ven = CoulombPotential::energies(samples_[id].sample_,atoms_);
        Eigen::VectorXd Ee = OneParticleEnergies::oneElectronEnergies(Te, Vee, Ven, Vnn_);
        TeStats_.add(Te,1);
        VeeStats_.add(Vee,1);
        VenStats_.add(Ven,1);
        EeStats_.add(Ee,1);
    }
    return count;
}

void EnergyCalculator::calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){
    using namespace YAML;
    using namespace Logger;

    yamlDocument_ << Key << "En" << Comment("[Eh]") << Value << EnStats_
                  << Key << "Clusters" << BeginSeq;

    size_t totalCount = 0;
    for (auto& cluster : clusteredGloballySimilarMaxima) {
        valueStats_.reset();
        SeeStats_.reset();
        TeStats_.reset();
        EeStats_.reset();
        VeeStats_.reset();
        VenStats_.reset();

        auto types = cluster[0].representativeReference().maximum().typesVector().asEigenVector();

        std::vector<ElectronsVector> structures;
        for (auto &simRefVector : cluster) {

            // Iterate over references being similar to the representative reference.

            for (const auto &ref : simRefVector.similarReferencesIterators()){
                totalCount += addReference(*ref);
            }
            structures.push_back(simRefVector.representativeReference().maximum());
        }
        printCluster(structures);
    }
    Logger::console->info("overall count {}", totalCount);
    assert(totalCount == samples_.size() && "The total count must match the sample size.");

    yamlDocument_ << EndSeq;
    assert(yamlDocument_.good());
}


// selects nWanted structures and prints the statistic data
void EnergyCalculator::printCluster(std::vector<ElectronsVector>& structures){

    size_t nWanted = 16;
    std::vector<ElectronsVector> selectedStructures;

    if(structures.size() < nWanted){
        for (auto& i : structures) selectedStructures.push_back(i);
    } else {
        size_t skipLength = (structures.size()-1) / (nWanted-1);

        for (size_t i = 0; i < nWanted; ++i)
            selectedStructures.push_back(structures[i*skipLength]);
    }

    yamlDocument_ << ClusterData(TeStats_.getTotalWeight(), selectedStructures, valueStats_, TeStats_, EeStats_,
                                 SeeStats_, VeeStats_, VenStats_);
}

YAML::Node EnergyCalculator::getYamlNode(){
    return YAML::Load(yamlDocument_.c_str());
}

std::string EnergyCalculator::getYamlDocumentString(){
    return std::string(yamlDocument_.c_str());
}
