// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <EnergyPartitioning.h>
#include "GraphAnalysis.h"
#include <Motifs.h>
#include <CoulombPotential.h>
#include <yaml-cpp/yaml.h>

using namespace testing;

class AEnergyPartitioningTest : public ::testing::Test {
public:

    Eigen::VectorXd Te;
    Eigen::MatrixXd Vee;
    Eigen::MatrixXd Ven;
    Eigen::MatrixXd Vnn;

    double totalEnergy;

    void SetUp() override {
        Te = Eigen::VectorXd::Constant(4,0.01);
        Vee = Eigen::MatrixXd::Constant(4,4,1);
        Ven = Eigen::MatrixXd::Constant(4,2,100);
        Vnn = Eigen::MatrixXd::Constant(2,2,10000);

        auto ne = Ven.rows();
        auto na = Ven.cols();

        totalEnergy = 0;
        for (int i = 0; i < ne; ++i) {
            totalEnergy += Te(i);

            for (int j = i+1; j < ne; ++j)
                totalEnergy += Vee(i,j);

            for (int k = 0; k < na; ++k)
                totalEnergy += Ven(i,k);
        }

        for (int k = 0; k < na; ++k) {
            for (int l = k+1; l < na; ++l)
                totalEnergy += Vnn(k,l);
        }

    }
};

TEST_F(AEnergyPartitioningTest, MolecularSelections) {
    auto selA = MolecularSelection({{0,1}},{{0}});
    auto selB = MolecularSelection({{2,3}},{{1}});

    auto res = EnergyPartitioning::MolecularSelectionBased::calculateInteractionEnergies(
            {selA, selB}, Te, Vee, Ven, Vnn);

    ASSERT_EQ(res.first.E.size(), 2);
    ASSERT_EQ(res.second.E.size(), 4);
    ASSERT_EQ(res.first.E[0], 2*0.01 + 1*1 + 2*100);
    ASSERT_EQ(res.first.E[1], 2*0.01 + 1*1 + 2*100);

    ASSERT_EQ(res.second.E(0,1),4*1 + 4*100 + 1*10000);
    ASSERT_EQ(res.second.E(0,0), 0);

    // matrix should be symmetric
    ASSERT_EQ(res.second.E(0,1), res.second.E(1,0));
    ASSERT_EQ(res.second.E(0,0), res.second.E(1,1));
}


TEST_F(AEnergyPartitioningTest, TwoPairs){

    MolecularGeometry mol = {
            AtomsVector({{Element::Li, {0, 0, 1}},
                         {Element::Li, {0, 0,-1}}}),
            ElectronsVector({{Spin::alpha, {0, 0, 1}},
                             {Spin::beta,  {0, 0, 1}},
                             {Spin::alpha, {0, 0,-1}},
                             {Spin::beta,  {0, 0,-1}},
                            })};

    Eigen::MatrixXb A(4, 4);
    //all pairs // Two pairs: (0-1,) (2-3)
    A <<
    0, 1, 1, 1, \
    1, 0, 1, 1, \
    1, 1, 0, 1, \
    1, 1, 1, 0;

    Motifs motifs(A, mol);

    ASSERT_THAT(motifs.motifs_[0].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(motifs.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifs_[1].electrons_.indices(), ElementsAre(2, 3));
    ASSERT_THAT(motifs.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifs_[1].type(), MotifType::Core);

    auto motifEnergies = EnergyPartitioning::MolecularSelectionBased::calculateInteractionEnergies(motifs, Te, Vee, Ven, Vnn);

    Eigen::VectorXd intraExpected = Eigen::VectorXd::Zero(motifs.motifs_.size());
    intraExpected << 201.02, 201.02;

    Eigen::MatrixXd interExpected = Eigen::MatrixXd::Zero(motifs.motifs_.size(), motifs.motifs_.size());
    interExpected << 0, 10404, 10404, 0;

    ASSERT_TRUE(motifEnergies.first.E.isApprox(intraExpected));
    ASSERT_TRUE(motifEnergies.second.E.isApprox(interExpected));

    ASSERT_EQ(motifEnergies.first.E.sum()+motifEnergies.second.E.sum()/2, totalEnergy);
};

TEST_F(AEnergyPartitioningTest, HydrogenMotif){
    MolecularGeometry mol = {
            AtomsVector({{Element::Li, {0, 0, 1}},
                         {Element::H, {0, 0,-1}}}),
            ElectronsVector({{Spin::alpha, {0, 0, 1}},
                             {Spin::beta,  {0, 0, 1}},
                             {Spin::alpha, {0, 0,-1}},
                             {Spin::beta,  {0, 0,-1.1}},
                            })};

    Eigen::MatrixXb A(4, 4);
    A << 0, 1, 0, 0, \
         1, 0, 0, 0, \
         0, 0, 0, 1, \
         0, 0, 1, 0;

    Motifs motifs(A, mol);

    ASSERT_THAT(motifs.motifs_[0].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(motifs.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifs_[1].electrons_.indices(), ElementsAre());
    ASSERT_THAT(motifs.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifs_[1].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifs_[2].electrons_.indices(), ElementsAre(2, 3));
    ASSERT_THAT(motifs.motifs_[2].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(motifs.motifs_[2].type(), MotifType::Valence);

    auto motifEnergies = EnergyPartitioning::MolecularSelectionBased::calculateInteractionEnergies(motifs, Te, Vee, Ven, Vnn);

    Eigen::VectorXd intraExpected = Eigen::VectorXd::Zero(motifs.motifs_.size());
    intraExpected << 201.02, 0.00, 1.02;

    Eigen::MatrixXd interExpected = Eigen::MatrixXd::Zero(motifs.motifs_.size(), motifs.motifs_.size());
    interExpected <<
    0, 10200, 204,
    10200, 0, 200,
    204, 200,0;

    ASSERT_TRUE(motifEnergies.first.E.isApprox(intraExpected));
    ASSERT_TRUE(motifEnergies.second.E.isApprox(interExpected));

    ASSERT_EQ(motifEnergies.first.E.sum()+motifEnergies.second.E.sum()/2, totalEnergy);
};

