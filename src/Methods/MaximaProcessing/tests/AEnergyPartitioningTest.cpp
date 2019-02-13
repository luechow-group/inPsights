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


TEST_F(AEnergyPartitioningTest, TwoPairs){

    MolecularGeometry mol = {
            AtomsVector({{Element::H, {0, 0, 1}},
                         {Element::H, {0, 0,-1}}}),
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

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre(2,3));
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, Te, Vee, Ven, Vnn);

    Eigen::VectorXd intraExpected = Eigen::VectorXd::Zero(motifs.motifVector_.size());
    intraExpected << 201.02, 201.02;

    Eigen::MatrixXd interExpected = Eigen::MatrixXd::Zero(motifs.motifVector_.size(),motifs.motifVector_.size());
    interExpected << 0, 10404, 0, 0;

    ASSERT_TRUE(motifEnergies.first.isApprox(intraExpected));
    ASSERT_TRUE(motifEnergies.second.isApprox(interExpected));

    ASSERT_EQ(motifEnergies.first.sum()+motifEnergies.second.sum(), totalEnergy);
};

TEST_F(AEnergyPartitioningTest, TripleAndIsolated){
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

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre(2));
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[2].electronIndices(), ElementsAre(3));
    ASSERT_THAT(motifs.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(motifs.motifVector_[2].type(), MotifType::Valence);

    auto motifEnergies = EnergyPartitioning::MotifBased::calculateInterationEnergies(motifs, Te, Vee, Ven, Vnn);

    Eigen::VectorXd intraExpected = Eigen::VectorXd::Zero(motifs.motifVector_.size());
    intraExpected << 201.02, 100.01, 0.01;

    Eigen::MatrixXd interExpected = Eigen::MatrixXd::Zero(motifs.motifVector_.size(),motifs.motifVector_.size());
    interExpected <<
    0, 10302, 102,
    0, 0, 101,
    0,0,0;

    ASSERT_TRUE(motifEnergies.first.isApprox(intraExpected));
    ASSERT_TRUE(motifEnergies.second.isApprox(interExpected));

    ASSERT_EQ(motifEnergies.first.sum()+motifEnergies.second.sum(), totalEnergy);
};

