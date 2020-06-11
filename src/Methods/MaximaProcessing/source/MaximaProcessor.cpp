/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <MaximaProcessor.h>
#include <MaximaProcessingSettings.h>
#include <ClusterData.h>
#include <MolecularGeometry.h>
#include <SpinCorrelation.h>
#include <CoulombPotential.h>
#include <EnergyPartitioning.h>
#include <VoxelCubeGeneration.h>
#include <VoxelCubeOverlapCalculation.h>
#include <spdlog/spdlog.h>
#include <SpinCorrelationValueHistogram.h>
#include <LocalParticleEnergiesCalculation.h>

MotifEnergyCalculator::Result MotifEnergyCalculator::partition(
        const Group &group,
        const std::vector<Sample> &samples,
        const Motifs& motifs) {

    VectorStatistics intraMotifEnergyStats;
    TriangularMatrixStatistics interMotifEnergyStats;

    partitionLowerLevels(group, samples, motifs, intraMotifEnergyStats, interMotifEnergyStats);

    return {intraMotifEnergyStats, interMotifEnergyStats};
}

void MotifEnergyCalculator::partitionLowerLevels(
        const Group &group,
        const std::vector<Sample> &samples,
        const Motifs& motifs,
        VectorStatistics& intraMotifEnergyStats,
        TriangularMatrixStatistics& interMotifEnergyStats) {

    if (group.isLeaf())
        partitionLowestLevel(group, samples, motifs, intraMotifEnergyStats, interMotifEnergyStats);
    else
        for (const auto &subgroup : group)
            partitionLowerLevels(subgroup, samples, motifs, intraMotifEnergyStats, interMotifEnergyStats);

}


void MotifEnergyCalculator::partitionLowestLevel(
        const Group &group,
        const std::vector<Sample> &samples,
        const Motifs& motifs,
        VectorStatistics& intraMotifEnergyStats,
        TriangularMatrixStatistics& interMotifEnergyStats) {

    auto finalNuclearPerm = group.representative()->nuclearPermutation();
    auto permutedNuclei =  group.representative()->nuclei();
    permutedNuclei.permute(finalNuclearPerm);

    for (auto id : group.representative()->sampleIds()) {
        Eigen::VectorXd Te = samples[id].kineticEnergies_;
        Eigen::MatrixXd Vee = CoulombPotential::energies(samples[id].sample_);
        Eigen::MatrixXd Ven = CoulombPotential::energies(samples[id].sample_, permutedNuclei);

        Eigen::MatrixXd Vnn = CoulombPotential::energies(permutedNuclei);

        auto motifEnergies = EnergyPartitioning::MotifBased::calculateInteractionEnergies(motifs, Te, Vee, Ven, Vnn);

        intraMotifEnergyStats.add(motifEnergies.first, 1);
        interMotifEnergyStats.add(motifEnergies.second, 1);
    }
}

MaximaProcessor::MaximaProcessor(YAML::Emitter& yamlDocument, const std::vector<Sample>& samples, const AtomsVector& atoms)
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
    valueStats_.add(value, count);
    SeeStats_.add(spinCorrelations_, count);

    auto permutedNuclei = reference.nuclei();
    permutedNuclei.permute(reference.nuclearPermutation());

    // Sample related statistics
    for (auto & id : reference.sampleIds()) {
        auto & electrons = samples_[id].sample_;
        Eigen::VectorXd Te = samples_[id].kineticEnergies_;
        Eigen::MatrixXd Vee = CoulombPotential::energies(electrons);
        Eigen::MatrixXd Ven = CoulombPotential::energies(electrons, permutedNuclei);

        Eigen::VectorXd Ee = EnergyPartitioning::ParticleBased::oneElectronEnergies(Te, Vee, Ven, Vnn_); // TODO wrong and deprecated -> remove
        Eigen::MatrixXd Ree = Metrics::positionalDistances(electrons.positionsVector());
        Eigen::MatrixXd Ren = Metrics::positionalDistances(electrons.positionsVector(), permutedNuclei.positionsVector());

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


// TODO every group should know its total count at merging already => add count as a member
size_t  MaximaProcessor::addAllReferences(const Group &group) {
    unsigned long totalCount = 0;
    if(group.isLeaf())
        totalCount += addReference(*group.representative());
    else
        for(const auto &subgroup : group)
            totalCount += addAllReferences(subgroup);

    return totalCount;
}

std::vector<ElectronsVector> MaximaProcessor::getAllRepresentativeMaxima(const Group &group) {
    std::vector<ElectronsVector> representativeMaxima;
    if(group.empty())
        representativeMaxima.emplace_back(group.representative()->maximum());
    else
        for(auto & i : group)
            representativeMaxima.emplace_back(i.representative()->maximum());

    return representativeMaxima;
}

void MaximaProcessor::calculateStatistics(const Group &maxima,
                                          const std::vector<std::vector<std::vector<size_t>>> & nucleiMergeLists,
                                          const std::vector<size_t> &nucleiIndices
                                          ){
    using namespace YAML;

    yamlDocument_ << Key << "Vnn" << Comment("[Eh]") << Value << VnnStats_
                  << Key << "En" << Comment("[Eh]") << Value << EnStats_
                  << Key << "Clusters" << BeginSeq;

    size_t totalCount = 0;
    double totalWeight = 0.0;

    SpinCorrelationValueHistogram spinCorrelationDistribution(12); // => 25 bns in total

    LocalParticleEnergiesCalculator localParticleEnergiesCalculator(
            samples_, atoms_, nucleiIndices,
            ParticleSelection::settings.maximalCount());

    for (auto& group : maxima) {

        valueStats_.reset();
        SeeStats_.reset();
        TeStats_.reset();
        EeStats_.reset();
        VeeStats_.reset();
        VenStats_.reset();
        EtotalStats_.reset();
        ReeStats_.reset();
        RenStats_.reset();


        totalCount += addAllReferences(group); // this sets all statistic objects internally
        auto structures = getAllRepresentativeMaxima(group);

        // SpinCorrelationDistribution
        spinCorrelationDistribution.addSpinStatistic(SeeStats_);
        
        // Motif analysis (requires spin correlation data)

        auto adjacencyMatrix = GraphAnalysis::filter(SeeStats_.mean().cwiseAbs(), MaximaProcessing::settings.motifThreshold());

        auto mol = MolecularGeometry(atoms_, group.representative()->maximum());
        Motifs motifs(adjacencyMatrix, mol);
        
        // merge motifs
        for(const auto& nucleiMergeList : nucleiMergeLists){
            auto motifMergeIndices = motifs.findMotifMergeIndices(mol, nucleiMergeList);

            motifs.mergeMotifs(motifMergeIndices);
        }

        // Motif energies
        auto [intraMotifEnergyStats, interMotifEnergyStats] = MotifEnergyCalculator::partition(group, samples_, motifs);

        // SEDs
        std::vector<VoxelCube> voxelCubes;
        if(VoxelCubeGeneration::settings.generateVoxelCubesQ())
            voxelCubes = VoxelCubeGeneration::fromCluster(group, samples_);

        // SED overlaps
        Eigen::MatrixXd overlaps;
        if(VoxelCubeOverlapCalculation::settings.calculateOverlapQ())
            overlaps = VoxelCubeOverlapCalculation::fromCluster(group, samples_);

        auto weight = double(TeStats_.getTotalWeight())/double(samples_.size());
        if(weight >= MaximaProcessing::settings.minimalClusterWeight.get()) {

            localParticleEnergiesCalculator.add(group);

            totalWeight += weight;

            ElectronsVector sampleAverage = {group.electronsVectorFromAveragedPositionsVector(group.averagedSamplePositionsVector(samples_))};

            yamlDocument_ << ClusterData(TeStats_.getTotalWeight(), structures, sampleAverage,
                    valueStats_, TeStats_, EeStats_,
                                         SeeStats_, VeeStats_, VenStats_,
                                         motifs, EtotalStats_, intraMotifEnergyStats, interMotifEnergyStats,
                                         ReeStats_, RenStats_, voxelCubes, overlaps);
        }
    }
    spdlog::info("Overall count {}. Some structures might be lost "
                 "by being classified as noise during density-based clustering.", totalCount);
    spdlog::info("Considered weight: {}, discarded weight: {}", totalWeight, 1.0 - totalWeight);

    yamlDocument_ << EndSeq;

    // add spinCorrelationDistribution
    auto hist = spinCorrelationDistribution.getHistogramVector();
    hist /= hist.sum();
    yamlDocument_ << Key << "SpinCorrelationDistribution" << Value << hist;

    yamlDocument_ << Key << "LocalParticleEnergiesCalculation"  << Value << localParticleEnergiesCalculator;



    assert(yamlDocument_.good());
}

YAML::Node MaximaProcessor::getYamlNode(){
    return YAML::Load(yamlDocument_.c_str());
}

std::string MaximaProcessor::getYamlDocumentString(){
    return std::string(yamlDocument_.c_str());
}
