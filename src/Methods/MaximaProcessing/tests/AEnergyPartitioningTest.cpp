//
// Created by Michael Heuer on 2019-02-03.
//

#include <gmock/gmock.h>
#include <EnergyPartitioning.h>
#include <GraphAnalysis.h>
#include <Motifs.h>
#include <CoulombPotential.h>
#include <yaml-cpp/yaml.h>

using namespace testing;

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

    Eigen::MatrixXd Vnn = Eigen::MatrixXd::Constant(4,4,10000);

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

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre(2,3));
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({0})], 201.02);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({0,1})], 10404);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({1})], 201.02);

    auto node = YAML::convert<MotifEnergies>::encode(motifEnergies);

    ASSERT_THAT(node[0]["MotifIds"].as<std::vector<long>>(), ElementsAre(0));
    ASSERT_THAT(node[1]["MotifIds"].as<std::vector<long>>(), ElementsAre(0,1));
    ASSERT_THAT(node[2]["MotifIds"].as<std::vector<long>>(), ElementsAre(1));

    MotifEnergies decodedMotifEnergies;
    YAML::convert<MotifEnergies>::decode(node, decodedMotifEnergies);

    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({0})], 201.02);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({0,1})], 10404);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({1})], 201.02);
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

    Motifs motifs(A, mol);

    SingleParticlesStatistics TeStats;
    IntraParticlesStatistics VeeStats;
    InterParticlesStatistics VenStats;

    TeStats.add(Eigen::VectorXd::Constant(mol.electrons().numberOfEntities(), 0.01));
    VeeStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.electrons().numberOfEntities(), 1));
    VenStats.add(Eigen::MatrixXd::Constant(mol.electrons().numberOfEntities(), mol.atoms().numberOfEntities(),100));

    Eigen::MatrixXd Vnn = Eigen::MatrixXd::Constant(4,4,10000);

    EnergyStatistics::ElectronicEnergy electronicEnergy(TeStats, VeeStats, VenStats);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, electronicEnergy, Vnn);

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre(2));
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[2].electronIndices(), ElementsAre(3));
    ASSERT_THAT(motifs.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(motifs.motifVector_[2].type(), MotifType::Valence);

    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({0})], 201.02);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({0,1})], 10302);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({0,2})], 102);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({1})], 100.01);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({1,2})], 101);
    ASSERT_EQ(motifEnergies.interactionEnergies()[std::vector<long>({2})], 0.01);


    auto node = YAML::convert<MotifEnergies>::encode(motifEnergies);

    ASSERT_THAT(node[0]["MotifIds"].as<std::vector<long>>(), ElementsAre(0));
    ASSERT_THAT(node[1]["MotifIds"].as<std::vector<long>>(), ElementsAre(0,1));
    ASSERT_THAT(node[2]["MotifIds"].as<std::vector<long>>(), ElementsAre(0,2));
    ASSERT_THAT(node[3]["MotifIds"].as<std::vector<long>>(), ElementsAre(1));
    ASSERT_THAT(node[4]["MotifIds"].as<std::vector<long>>(), ElementsAre(1,2));
    ASSERT_THAT(node[5]["MotifIds"].as<std::vector<long>>(), ElementsAre(2));


    MotifEnergies decodedMotifEnergies;
    YAML::convert<MotifEnergies>::decode(node, decodedMotifEnergies);

    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({0})], 201.02);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({0,1})], 10302);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({0,2})], 102);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({1})], 100.01);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({1,2})], 101);
    ASSERT_EQ(decodedMotifEnergies.interactionEnergies()[std::vector<long>({2})], 0.01);
};

