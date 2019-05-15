//
// Created by heuer on 12.12.18.
//
#include <MaximaProcessor.h>
#include <MaximaProcessingSettings.h>
#include <ClusterData.h>
#include <MolecularGeometry.h>
#include <SpinCorrelation.h>
#include <CoulombPotential.h>
#include <EnergyPartitioning.h>
#include <VoxelCubeGeneration.h>
#include <EnergyPartitioning.h>
#include <spdlog/spdlog.h>

MaximaProcessor::MaximaProcessor(YAML::Emitter& yamlDocument, const std::vector<Sample>& samples, AtomsVector atoms)
        :
        yamlDocument_(yamlDocument),
        samples_(samples),
        atoms_(atoms),
        Vnn_(CoulombPotential::energies<Element>(atoms_))
{
    VnnStats_.add(Vnn_);
    EnStats_.add(EnergyPartitioning::ParticleBased::oneAtomEnergies(Vnn_));
}

unsigned long MaximaProcessor::addReference(const Reference &reference) {
    auto count = unsigned(reference.count());

    Eigen::VectorXd value = Eigen::VectorXd::Constant(1, reference.value());
    Eigen::MatrixXd spinCorrelations_ = SpinCorrelation::spinCorrelations(reference.maximum().typesVector()).cast<double>();
    // Maximum related statistics
    std::cout << reference.maximum().typesVector() << std::endl;
    std::cout << spinCorrelations_ << std::endl;
    valueStats_.add(value, count);
    SeeStats_.add(spinCorrelations_, count);

    // Sample related statistics
    for (auto & id : reference.sampleIds()) {
        Eigen::VectorXd Te = samples_[id].kineticEnergies_;
        Eigen::MatrixXd Vee = CoulombPotential::energies(samples_[id].sample_);
        Eigen::MatrixXd Ven = CoulombPotential::energies(samples_[id].sample_,atoms_);
        Eigen::VectorXd Ee = EnergyPartitioning::ParticleBased::oneElectronEnergies(Te, Vee, Ven, Vnn_);
        Eigen::MatrixXd Ree = Metrics::positionalDistances(samples_[id].sample_.positionsVector());
        Eigen::MatrixXd Ren = Metrics::positionalDistances(samples_[id].sample_.positionsVector(), atoms_.positionsVector());

        TeStats_.add(Te,1);
        VeeStats_.add(Vee,1);
        VenStats_.add(Ven,1);
        EeStats_.add(Ee,1);
        ReeStats_.add(Ree,1);
        RenStats_.add(Ren,1);

        Eigen::VectorXd Etot(1);
        Etot << EnergyPartitioning::calculateTotalEnergy(Te,Vee,Ven, Vnn_);
        EtotalStats_.add(Etot);
    }
    return count;
}

void MaximaProcessor::doMotifBasedEnergyPartitioning(const Group& group) {

    auto allSampleIds = group.allSampleIds();
    
    // Sample related statistics
    for (auto id : allSampleIds) {
        Eigen::VectorXd Te = samples_[id].kineticEnergies_;
        Eigen::MatrixXd Vee = CoulombPotential::energies(samples_[id].sample_);
        Eigen::MatrixXd Ven = CoulombPotential::energies(samples_[id].sample_,atoms_);

        auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs_, Te, Vee, Ven, Vnn_);

        intraMotifEnergyStats_.add(motifEnergies.first,1);
        interMotifEnergyStats_.add(motifEnergies.second,1);
    }
}

size_t  MaximaProcessor::addAllReferences(const Group &group) {
    std::cout << group << std::endl;
    unsigned long totalCount = 0;
    if(group.isLeaf())
        totalCount += addReference(*group.representative());
    else
        for(const auto &subgroup : group)
            totalCount += addAllReferences(subgroup);

    return totalCount;
}

std::vector<ElectronsVector> MaximaProcessor::getAllRepresentativeMaxima(const Group &group) {
    if(group.isLeaf()) {
        return {}; // leaves are not printed directly (only as representative structures one layer above)
    } else {
        if(group.front().isLeaf()) {
            return {group.representative()->maximum()};
        } else {
            std::vector<ElectronsVector> representativeMaxima;
            for(auto & i : group) {
                auto res = getAllRepresentativeMaxima(i);
                representativeMaxima.insert(representativeMaxima.end(), res.begin(), res.end());
            }
            return representativeMaxima;
        }
    }
}

void MaximaProcessor::calculateStatistics(const Group &maxima){
    using namespace YAML;

    yamlDocument_ << Key << "En" << Comment("[Eh]") << Value << EnStats_
                  << Key << "Clusters" << BeginSeq;

    size_t totalCount = 0;
    for (auto& group : maxima) {

        valueStats_.reset();
        SeeStats_.reset();
        TeStats_.reset();
        EeStats_.reset();
        VeeStats_.reset();
        VenStats_.reset();
        EtotalStats_.reset();
        intraMotifEnergyStats_.reset();
        interMotifEnergyStats_.reset();
        ReeStats_.reset();
        RenStats_.reset();


        totalCount += addAllReferences(group);
        auto structures = getAllRepresentativeMaxima(group);
        
        // Motif analysis (requires spin correlation data)
        auto adjacencyMatrix = GraphAnalysis::filter(SeeStats_.mean().cwiseAbs(), 1.00);
        motifs_ = Motifs(adjacencyMatrix, MolecularGeometry(atoms_, group.representative()->maximum()));
        
        doMotifBasedEnergyPartitioning(group);

        std::vector<VoxelCube> voxelCubes;
        if(VoxelCubeGeneration::settings.generateVoxelCubesQ())
            voxelCubes = VoxelCubeGeneration::fromCluster(group, samples_);


        if(TeStats_.getTotalWeight() > MaximaProcessing::settings.minimalClusterSize.get()) {
            yamlDocument_ << ClusterData(TeStats_.getTotalWeight(), structures, valueStats_, TeStats_, EeStats_,
                                         SeeStats_, VeeStats_, VenStats_,
                                         motifs_, EtotalStats_, intraMotifEnergyStats_, interMotifEnergyStats_,
                                         ReeStats_, RenStats_, voxelCubes);
        }
    }
    spdlog::info("overall count {}", totalCount);
    assert(totalCount == samples_.size() && "The total count must match the sample size.");

    yamlDocument_ << EndSeq;
    assert(yamlDocument_.good());
}

YAML::Node MaximaProcessor::getYamlNode(){
    return YAML::Load(yamlDocument_.c_str());
}

std::string MaximaProcessor::getYamlDocumentString(){
    return std::string(yamlDocument_.c_str());
}
