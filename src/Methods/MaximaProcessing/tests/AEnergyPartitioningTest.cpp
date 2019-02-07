//
// Created by Michael Heuer on 2019-02-03.
//

#include <gmock/gmock.h>
#include <EnergyPartitioning.h>
#include <GraphAnalysis.h>
#include <Motifs.h>
#include <CoulombPotential.h>
#include <yaml-cpp/yaml.h>

TEST(AEnergyPartitioningTest, TwoPairs){


    MolecularGeometry mol = {
            AtomsVector({{Element::H, {0, 0, 1}},
                         {Element::H, {0, 0,-1}}}),
            ElectronsVector({{Spin::alpha, {0, 0, 1}},
                             {Spin::beta,  {0, 0, 1}},
                             {Spin::alpha, {0, 0,-1}},
                             {Spin::beta,  {0, 0,-1}},
                            })};

    SingleParticlesStatistics TeStats;
    IntraParticlesStatistics VeeStats;
    InterParticlesStatistics VenStats;

    TeStats.add(Eigen::VectorXd::Constant(mol.electrons().numberOfEntities(), 0.01));
    VeeStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.electrons().numberOfEntities(), 1));
    VenStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.atoms().numberOfEntities(),100));

    Eigen::MatrixXd Vnn(4,4);
    Vnn <<
    10000,10000,10000,10000,\
    10000,10000,10000,10000,\
    10000,10000,10000,10000,\
    10000,10000,10000,10000;

    Eigen::MatrixXb A(4, 4);
    //all pairs // Two pairs: (0-1,) (2-3)
    A <<
    0, 1, 1, 1, \
    1, 0, 1, 1, \
    1, 1, 0, 1, \
    1, 1, 1, 0;

    Motifs motifs(A, mol);
    EnergyStatistics::ElectronicEnergy electronicEnergy(TeStats, VeeStats, VenStats);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, electronicEnergy, Vnn);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({0,1}), Motif({0,1}))], 201.02);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({0,1}), Motif({2,3}))], 10404);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({2,3}), Motif({2,3}))], 201.02);
};



TEST(AEnergyPartitioningTest, TripleAndIsolated){
    MolecularGeometry mol = {
            AtomsVector({{Element::H, {0, 0, 1}},
                         {Element::H, {0, 0,-1}}}),
            ElectronsVector({{Spin::alpha, {0, 0, 1}},
                             {Spin::beta,  {0, 0, 1}},
                             {Spin::alpha, {0, 0,-1}},
                             {Spin::beta,  {0, 0,-1.1}},
                            })};

    Eigen::MatrixXb A(4, 4);
    A << 0, 1, 1, 0, \
         1, 0, 1, 0, \
         1, 1, 0, 0, \
         0, 0, 0, 0;


    SingleParticlesStatistics TeStats;
    IntraParticlesStatistics VeeStats;
    InterParticlesStatistics VenStats;

    TeStats.add(Eigen::VectorXd::Constant(mol.electrons().numberOfEntities(), 0.01));
    VeeStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.electrons().numberOfEntities(), 1));
    VenStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.atoms().numberOfEntities(),100));

    Eigen::MatrixXd Vnn(4,4);
    Vnn <<
    10000,10000,10000,10000,\
    10000,10000,10000,10000,\
    10000,10000,10000,10000,\
    10000,10000,10000,10000;

    Motifs motifs(A, mol);
    EnergyStatistics::ElectronicEnergy electronicEnergy(TeStats, VeeStats, VenStats);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, electronicEnergy, Vnn);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({0,1}), Motif({0,1}))], 201.02);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({0,1}), Motif({2}))], 10302);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({0,1}), Motif({3}))], 102);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({2}), Motif({2}))], 100.01);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({2}), Motif({3}))], 101);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::make_pair(Motif({3}), Motif({3}))], 0.01);
};