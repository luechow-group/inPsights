//
// Created by Michael Heuer on 2019-02-03.
//

#include <gmock/gmock.h>
#include <EnergyPartitioning.h>
#include <GraphAnalysis.h>
#include <MotifAnalysis.h>

TEST(AEnergyPartitioningTest, Test){
    auto electronsCount = 4;

    SingleParticlesStatistics TeStats;
    IntraParticlesStatistics VeeStats;
    InterParticlesStatistics VenStats;


    TeStats.add(Eigen::VectorXd::Constant(electronsCount,0.01));
    VeeStats.add(Eigen::MatrixXd::Constant(electronsCount,electronsCount,1));
    VenStats.add(Eigen::MatrixXd::Constant(electronsCount,electronsCount,100));

    Eigen::MatrixXb  A(4, 4);
    // Two pairs: (0-1,) (2-3)
    A <<
    0, 1, 0, 0, \
    1, 0, 0, 0, \
    0, 0, 0, 1, \
    0, 0, 1, 0;

    MotifAnalysis::Motifs motifs(A);

    EnergyStatistics::ElectronicEnergy electronicEnergy(TeStats, VeeStats, VenStats);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, electronicEnergy);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(
            MotifAnalysis::Motif({0,1}), MotifAnalysis::Motif({0,1}))], 401.01);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(
            MotifAnalysis::Motif({0,1}), MotifAnalysis::Motif({2,3}))], 4);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(
            MotifAnalysis::Motif({2,3}), MotifAnalysis::Motif({2,3}))], 401.01);

    /* TODO
     * KERNELEKTRONEN PER DISTANZKRITERIUM DEN ATOMEN ZUWEISEN
     * KERN-WW-ENERGIEN MIT VERRECHNEN
     */

};