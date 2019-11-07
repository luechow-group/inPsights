/* Copyright (C) 2019 Michael Heuer.
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

#include <gmock/gmock.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <limits>
#include <random>
#include <Metrics.h>
#include <Hungarian.h>
#include <Combinatorics.h>

using namespace testing;
using namespace SOAP;

class ABestMatchSimilarityTest : public ::testing::Test {
public:
    double distanceTolerance, soapThreshold, shakeSoapThreshold, eps;

    void SetUp() override {
        //spdlog::set_level(spdlog::level::off);

        distanceTolerance = 0.1;
        soapThreshold = 1.0;
        shakeSoapThreshold = 0.90;
        eps = 1E-9;

        Radial::settings.nmax = 2;
        Radial::settings.sigmaAtom = 2.0;
        Angular::settings.lmax = 2;
        Cutoff::settings.radius = 8.0;
        Cutoff::settings.width = 2.0;
        General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    };

    void routine(const MolecularGeometry & A, const MolecularGeometry & B,
                 const std::vector<Eigen::VectorXi> &expectedPermutationIndices,
                 double distTolerance, double soapThresh,bool lessThan = false){

        ParticleKit::create(A);

        auto specA = MolecularSpectrum(A);
        auto specB = MolecularSpectrum(B);

        auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distTolerance, soapThresh);

        for (auto& i : results)
            std::cout << i.metric << ", " << i.permutation.indices().transpose() << std::endl;
        for (auto& i : expectedPermutationIndices)
            std::cout << "   " << i.transpose() << std::endl;

            ASSERT_EQ(results.size(), expectedPermutationIndices.size());
        for (auto& i : results) {
            if(!lessThan)
                ASSERT_NEAR(i.metric, 1.0, eps);
            else
                ASSERT_LT(i.metric, 1.0);
            ASSERT_THAT(i.permutation.indices(), AnyOfArray(expectedPermutationIndices));
        }
    }
};


TEST_F(ABestMatchSimilarityTest, PrermuteEnvironmentsToLabSystem) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::BH3::ionicMinimal;
    auto A = TestMolecules::BH3::ionicMinimalRotatedPermuted;
    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto environmentalSimilarities = BestMatch::SOAPSimilarity::calculateEnvironmentalSimilarityMatrix(specA,specB);

    auto fromKitA = ParticleKit::fromKitPermutation(A.electrons());
    auto fromKitB = ParticleKit::fromKitPermutation(B.electrons());

    auto environmentalSimilaritiesInLabSystem = fromKitA*environmentalSimilarities*fromKitB;

    auto eps = std::numeric_limits<double>::epsilon()*10;

    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(0,0), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(0,1), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(1,2), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(1,3), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(2,0), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(2,1), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(3,2), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(3,3), 1.0, eps);
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndices) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,0,1,0,\
    0,1,0,1,\
    1,0,1,0,\
    0,1,0,1;

    auto bestMatch = Hungarian<double>::findMatching(mat, Matchtype::MAX);

    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(mat.rows());
    expectedPerm.setIdentity();
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto    result = BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(mat,bestMatch,1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_EQ(result[0][0].first, 0);
    ASSERT_EQ(result[0][1].first, 2);

    ASSERT_EQ(result[1][0].first, 1);
    ASSERT_EQ(result[1][1].first, 3);
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndices2) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,1,0,0,\
    1,1,0,0,\
    0,0,1,1,\
    0,0,1,1;

    auto bestMatch = Hungarian<double>::findMatching(mat, Matchtype::MAX);

    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(mat.rows());
    expectedPerm.setIdentity();
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto result = BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(mat,bestMatch,1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_EQ(result[0][0].first, 0);
    ASSERT_EQ(result[0][1].first, 1);

    ASSERT_EQ(result[1][0].first, 2);
    ASSERT_EQ(result[1][1].first, 3);
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndices3) {
    Eigen::MatrixXd mat(6,6);
    mat <<
    1,1,0,0,0,0,\
    0,0,1,0,0,0,\
    1,1,0,0,0,0,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1;

    auto bestMatch = Hungarian<double>::findMatching(mat, Matchtype::MAX);

    Eigen::VectorXi expectedPermIndices(mat.rows());
    expectedPermIndices << 0,2,1,3,4,5;
    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(expectedPermIndices);
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto result = BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(mat,bestMatch,1.0);

    ASSERT_EQ(result.size(),3);

    ASSERT_EQ(result[0].size(),2);
    ASSERT_EQ(result[0][0].first, 0);
    ASSERT_EQ(result[0][0].second,0);
    ASSERT_EQ(result[0][1].first, 2);
    ASSERT_EQ(result[0][1].second,1);

    ASSERT_EQ(result[1].size(),1);
    ASSERT_EQ(result[1][0].first, 1);
    ASSERT_EQ(result[1][0].second,2);

    ASSERT_EQ(result[2].size(),3);
    ASSERT_EQ(result[2][0].first, 3);
    ASSERT_EQ(result[2][0].second,3);
    ASSERT_EQ(result[2][1].first, 4);
    ASSERT_EQ(result[2][1].second,4);
    ASSERT_EQ(result[2][2].first, 5);
    ASSERT_EQ(result[2][2].second,5);
}


TEST_F(ABestMatchSimilarityTest, DISABLED_VarySimilarEnvironments) {
    //TODO
}

TEST_F(ABestMatchSimilarityTest, H4_Chemical) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4_Chemical_Shaked) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    ParticleKit::create(A);
    
    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for(auto seed : std::vector<unsigned long>{0,randomSeed}) {
        auto rng = std::default_random_engine(seed);

        while (Metrics::positionalNormsVectorNorm<Eigen::Dynamic, 2>(
                A.electrons().positionsVector(),
                B.electrons().positionsVector()) == 0.0)
            A.electrons().positionsVector().shake(distanceTolerance / 2.0, rng);

        std::vector<Eigen::VectorXi> permsIndices(4, Eigen::VectorXi(B.electrons().numberOfEntities()));
        permsIndices[0] << 0, 1, 2, 3;
        permsIndices[1] << 0, 1, 3, 2;
        permsIndices[2] << 1, 0, 2, 3;
        permsIndices[3] << 1, 0, 3, 2;

        routine(A, B, permsIndices, distanceTolerance, shakeSoapThreshold, true);
    }
}

TEST_F(ABestMatchSimilarityTest, H4_Alchemical) {
    General::settings.mode = General::Mode::alchemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4_Alchemical_Shaked) {
    General::settings.mode = General::Mode::alchemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    
    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for(auto seed : std::vector<unsigned long>{0,randomSeed}) {
        auto rng = std::default_random_engine(seed);

        while (Metrics::positionalNormsVectorNorm<Eigen::Dynamic, 2>(
                A.electrons().positionsVector(),
                B.electrons().positionsVector()) == 0.0)
            A.electrons().positionsVector().shake(distanceTolerance / 2.0, rng);

        std::vector<Eigen::VectorXi> permsIndices(4, Eigen::VectorXi(B.electrons().numberOfEntities()));
        permsIndices[0] << 0, 1, 2, 3;
        permsIndices[1] << 0, 1, 3, 2;
        permsIndices[2] << 1, 0, 2, 3;
        permsIndices[3] << 1, 0, 3, 2;

        routine(A, B, permsIndices, distanceTolerance, shakeSoapThreshold, true);
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_BoraneLEO1) {
    General::settings.mode = General::Mode::chemical;

    auto nuclei = TestMolecules::BH3::nuclei;
    const MolecularGeometry A = {
            TestMolecules::BH3::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                {Spin::alpha, nuclei.atoms().positionsVector()[3]}
            })};

    std::vector<Eigen::VectorXi> permsIndices(6, Eigen::VectorXi(A.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2;
    permsIndices[1] << 0,2,1;
    permsIndices[2] << 1,0,2;
    permsIndices[3] << 1,2,0;
    permsIndices[4] << 2,0,1;
    permsIndices[5] << 2,1,0;

    routine(A,A,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, BH3Covalent_Chemical) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    auto nuclei = BH3::nuclei;

    const MolecularGeometry A = {
            nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                {Spin::beta,inbetween(nuclei.atoms(),{0,1},0.25)},
                {Spin::beta,inbetween(nuclei.atoms(),{0,2},0.25)},
                {Spin::beta,inbetween(nuclei.atoms(),{0,3},0.25)}
            })};

    std::vector<Eigen::VectorXi> permsIndices(6, Eigen::VectorXi(A.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,2,1,3,5,4;
    permsIndices[2] << 1,0,2,4,3,5;
    permsIndices[3] << 1,2,0,4,5,3;
    permsIndices[4] << 2,0,1,5,3,4;
    permsIndices[5] << 2,1,0,5,4,3;

    routine(A,A,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, BH3IonicSelf_Chemical) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    auto nuclei = BH3::nuclei;

    const MolecularGeometry A = {
            nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                {Spin::beta,inbetween(nuclei.atoms(),{0,1},0.25)},
                {Spin::beta,inbetween(nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                {Spin::beta,inbetween(nuclei.atoms(),{0,3},0.25)}
            })};

    std::vector<Eigen::VectorXi> permsIndices(6, Eigen::VectorXi(A.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,4,2,5,1,3;
    permsIndices[2] << 1,0,3,2,4,5;
    permsIndices[3] << 1,4,3,5,0,2;
    permsIndices[4] << 4,0,5,2,1,3;
    permsIndices[5] << 4,1,5,3,0,2;

    routine(A,A,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, BH3_ThreeIndependentUnsimilarEnvironments_Chemical) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    auto nuclei = BH3::nuclei.atoms();

    const MolecularGeometry B = {
            nuclei,
            ElectronsVector({
                {Spin::alpha, nuclei.positionsVector()[1]},
                {Spin::alpha, nuclei.positionsVector()[2]},
                {Spin::beta,inbetween(nuclei,{0,1},0.25)}
                            })};

    const MolecularGeometry A = {
            nuclei,
            ElectronsVector({
                {Spin::alpha, nuclei.positionsVector()[1]},
                {Spin::alpha, nuclei.positionsVector()[2]},
                {Spin::beta,inbetween(nuclei,{0,2},0.25)}
            })};

    std::vector<Eigen::VectorXi> permsIndices(1, Eigen::VectorXi(A.electrons().numberOfEntities()));
    permsIndices[0] << 1,0,2;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}


TEST_F(ABestMatchSimilarityTest, BH3Ionic_Chemical) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    auto nuclei = BH3::nuclei.atoms();

    const MolecularGeometry B = {
            nuclei,
            ElectronsVector({
                {Spin::alpha, nuclei.positionsVector()[2]},
                {Spin::alpha, nuclei.positionsVector()[3]},
                 {Spin::beta , inbetween(nuclei,{0,2},0.25)},
                 {Spin::beta , inbetween(nuclei,{0,3},0.25)}
            })};

    const MolecularGeometry A = {
            nuclei,
            ElectronsVector({
                {Spin::alpha, nuclei.positionsVector()[1]},
                {Spin::alpha, nuclei.positionsVector()[2]},
                {Spin::beta , inbetween(nuclei,{0,1},0.25)},
                {Spin::beta , inbetween(nuclei,{0,2},0.25)}
            })};

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3;
    permsIndices[1] << 1,0, 3,2;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}

/*TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_EthaneIonic) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::Ethane::doublyIonicMinimal;
    auto A = TestMolecules::Ethane::doublyIonicMinimalDoublyRotated;

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3, 4,5, 6,7;
    permsIndices[1] << 1,0, 3,2, 5,4, 7,6;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}*/

TEST_F(ABestMatchSimilarityTest, EthaneIonic) {
    General::settings.mode = General::Mode::chemical;

    MolecularGeometry B = {
            TestMolecules::Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[2]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[2]},
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[5]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[5]}}
            )};
    MolecularGeometry A = {
            TestMolecules::Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[3]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[3]}}
            )};

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3;
    permsIndices[1] << 2,3, 0,1;

    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, EthaneSinglyIonicMinimal) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    MolecularGeometry B = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
                {Spin::alpha /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
                {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                {Spin::beta/*16*/,inbetween(Ethane::nuclei.atoms(),{1,6},0.25)},
                {Spin::beta/*17*/,Ethane::nuclei.atoms().positionsVector()[7]}
                            })};
    MolecularGeometry A = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                 {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                 {Spin::alpha /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
                 {Spin::alpha /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
                 {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                 {Spin::beta/*16*/,Ethane::nuclei.atoms().positionsVector()[6]},
                 {Spin::beta/*17*/,inbetween(Ethane::nuclei.atoms(),{1,7},0.25)}
            })};

    std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permIndices[0] << 1, 2, 0, 4, 5, 3; // 120° rotation around z
    permIndices[1] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, EthaneSinglyIonic) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    MolecularGeometry B = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
            // C Cores
            {Spin::alpha/* 0*/,Ethane::nuclei.atoms().positionsVector()[0]},
            {Spin::beta /* 1*/,Ethane::nuclei.atoms().positionsVector()[0]},
            {Spin::alpha/* 2*/,Ethane::nuclei.atoms().positionsVector()[1]},
            {Spin::beta /* 3*/,Ethane::nuclei.atoms().positionsVector()[1]},
            // H Cores
            {Spin::alpha/* 4*/,Ethane::nuclei.atoms().positionsVector()[2]},
            {Spin::alpha/* 5*/,Ethane::nuclei.atoms().positionsVector()[3]},
            {Spin::alpha/* 6*/,Ethane::nuclei.atoms().positionsVector()[4]},
            {Spin::beta /* 7*/,Ethane::nuclei.atoms().positionsVector()[5]},
            {Spin::beta /* 8*/,Ethane::nuclei.atoms().positionsVector()[6]},
            {Spin::beta /* 9*/,Ethane::nuclei.atoms().positionsVector()[7]},
            // C-C Bond
            {Spin::alpha/*10*/,inbetween(Ethane::nuclei.atoms(),{0,1},0.25)},
            {Spin::beta /*11*/,inbetween(Ethane::nuclei.atoms(),{1,0},0.25)},
            // C-H Bonds
            {Spin::beta /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
            {Spin::beta /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
            {Spin::beta /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
            {Spin::alpha/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
            {Spin::alpha/*16*/,inbetween(Ethane::nuclei.atoms(),{1,6},0.25)},
            {Spin::alpha/*17*/,Ethane::nuclei.atoms().positionsVector()[7]}
            })};
    MolecularGeometry A = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
            // C Cores
            {Spin::alpha/* 0*/,Ethane::nuclei.atoms().positionsVector()[0]},
            {Spin::beta /* 1*/,Ethane::nuclei.atoms().positionsVector()[0]},
            {Spin::alpha/* 2*/,Ethane::nuclei.atoms().positionsVector()[1]},
            {Spin::beta /* 3*/,Ethane::nuclei.atoms().positionsVector()[1]},
            // H Cores
            {Spin::alpha/* 4*/,Ethane::nuclei.atoms().positionsVector()[2]},
            {Spin::alpha/* 5*/,Ethane::nuclei.atoms().positionsVector()[3]},
            {Spin::alpha/* 6*/,Ethane::nuclei.atoms().positionsVector()[4]},
            {Spin::beta /* 7*/,Ethane::nuclei.atoms().positionsVector()[5]},
            {Spin::beta /* 8*/,Ethane::nuclei.atoms().positionsVector()[6]},
            {Spin::beta /* 9*/,Ethane::nuclei.atoms().positionsVector()[7]},
            // C-C Bond
            {Spin::alpha/*10*/,inbetween(Ethane::nuclei.atoms(),{0,1},0.25)},
            {Spin::beta /*11*/,inbetween(Ethane::nuclei.atoms(),{1,0},0.25)},
            // C-H Bonds
            {Spin::beta /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
            {Spin::beta /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
            {Spin::beta /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
            {Spin::alpha/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
            {Spin::alpha/*16*/,Ethane::nuclei.atoms().positionsVector()[6]},
            {Spin::alpha/*17*/,inbetween(Ethane::nuclei.atoms(),{1,7},0.25)}
            })};

    std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permIndices[0] << 0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15; // 120° rotation around z
    permIndices[1] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, DISABLED_EthaneSinglyIonicPermutedMinimal) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    MolecularGeometry B = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
                {Spin::alpha /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
                {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                {Spin::beta/*16*/,inbetween(Ethane::nuclei.atoms(),{1,6},0.25)},
                {Spin::beta/*17*/,Ethane::nuclei.atoms().positionsVector()[7]}
                            })};
    MolecularGeometry A = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
                {Spin::alpha /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
                {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                {Spin::beta/*16*/,Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::beta/*17*/,inbetween(Ethane::nuclei.atoms(),{1,7},0.25)}
            })};


    std::vector<int> indices(A.electrons().numberOfEntities());
    std::iota(indices.begin(), indices.end(), 0);

    Combinatorics::Permutations<int> indexSwapPerms(indices);

    for(auto indexSwap : indexSwapPerms) {
        //Eigen::VectorXi indexSwap(A.electrons().numberOfEntities());
        ////indexSwap << 0, 1, 3, 2, 4, 5; // works
        ////indexSwap << 0, 1, 4, 2, 3, 5; // works
        //indexSwap << 0, 1, 2, 4, 3, 5; // works
        ////indexSwap << 1, 0, 2, 3, 4, 5; // works
        ////indexSwap << 0, 1, 2, 3, 5, 4; // works
        //Eigen::PermutationMatrix<Eigen::Dynamic> indexSwapPerm(indexSwap);

        Eigen::Map<Eigen::VectorXi> v(indexSwap.data(),indexSwap.size());
        std::cout << "new perm to try:" << v.transpose() << std::endl;
        Eigen::PermutationMatrix<Eigen::Dynamic> indexSwapPerm(v);
        std::cout << "new perm to try:" << indexSwapPerm.indices().transpose() << std::endl;

        auto Acopy = A;
        Acopy.electrons().permute(indexSwapPerm);

        std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
        permIndices[0] << 1, 2, 0, 4, 5, 3; // 120° rotation around z
        permIndices[1] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane

        permIndices[0] = indexSwapPerm * permIndices[0];
        permIndices[1] = indexSwapPerm * permIndices[1];

        routine(Acopy, B, permIndices, distanceTolerance, soapThreshold);
    }
}


TEST_F(ABestMatchSimilarityTest, EthaneDoublyIonicAntiMinimal) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    MolecularGeometry B = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
                {Spin::alpha /*14*/,Ethane::nuclei.atoms().positionsVector()[4]},
                {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                {Spin::beta/*16*/,inbetween(Ethane::nuclei.atoms(),{1,6},0.25)},
                {Spin::beta/*17*/,Ethane::nuclei.atoms().positionsVector()[7]}
            })};
    MolecularGeometry A = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
                {Spin::alpha /*13*/,Ethane::nuclei.atoms().positionsVector()[3]},
                {Spin::alpha /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
                {Spin::beta/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
                {Spin::beta/*16*/,Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::beta/*17*/,inbetween(Ethane::nuclei.atoms(),{1,7},0.25)}
            })};

    std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permIndices[0] << 1, 2, 0, 4, 5, 3; // 120° rotation around z
    permIndices[1] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, EthaneDoublyIonicAnti) {
    General::settings.mode = General::Mode::chemical;

    using namespace TestMolecules;
    MolecularGeometry B = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
               // C Cores
               {Spin::alpha/* 0*/,Ethane::nuclei.atoms().positionsVector()[0]},
               {Spin::beta /* 1*/,Ethane::nuclei.atoms().positionsVector()[0]},
               {Spin::alpha/* 2*/,Ethane::nuclei.atoms().positionsVector()[1]},
               {Spin::beta /* 3*/,Ethane::nuclei.atoms().positionsVector()[1]},
               // H Cores
               {Spin::alpha/* 4*/,Ethane::nuclei.atoms().positionsVector()[2]},
               {Spin::alpha/* 5*/,Ethane::nuclei.atoms().positionsVector()[3]},
               {Spin::alpha/* 6*/,Ethane::nuclei.atoms().positionsVector()[4]},
               {Spin::beta /* 7*/,Ethane::nuclei.atoms().positionsVector()[5]},
               {Spin::beta /* 8*/,Ethane::nuclei.atoms().positionsVector()[6]},
               {Spin::beta /* 9*/,Ethane::nuclei.atoms().positionsVector()[7]},
               // C-C Bond
               {Spin::alpha/*10*/,inbetween(Ethane::nuclei.atoms(),{0,1},0.25)},
               {Spin::beta /*11*/,inbetween(Ethane::nuclei.atoms(),{1,0},0.25)},
               // C-H Bonds
               {Spin::beta /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
               {Spin::beta /*13*/,inbetween(Ethane::nuclei.atoms(),{0,3},0.25)},
               {Spin::beta /*14*/,Ethane::nuclei.atoms().positionsVector()[4]},
               {Spin::alpha/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
               {Spin::alpha/*16*/,inbetween(Ethane::nuclei.atoms(),{1,6},0.25)},
               {Spin::alpha/*17*/,Ethane::nuclei.atoms().positionsVector()[7]}
            })};
    MolecularGeometry A = {
            Ethane::nuclei.atoms(),
            ElectronsVector({
               // C Cores
               {Spin::alpha/* 0*/,Ethane::nuclei.atoms().positionsVector()[0]},
               {Spin::beta /* 1*/,Ethane::nuclei.atoms().positionsVector()[0]},
               {Spin::alpha/* 2*/,Ethane::nuclei.atoms().positionsVector()[1]},
               {Spin::beta /* 3*/,Ethane::nuclei.atoms().positionsVector()[1]},
               // H Cores
               {Spin::alpha/* 4*/,Ethane::nuclei.atoms().positionsVector()[2]},
               {Spin::alpha/* 5*/,Ethane::nuclei.atoms().positionsVector()[3]},
               {Spin::alpha/* 6*/,Ethane::nuclei.atoms().positionsVector()[4]},
               {Spin::beta /* 7*/,Ethane::nuclei.atoms().positionsVector()[5]},
               {Spin::beta /* 8*/,Ethane::nuclei.atoms().positionsVector()[6]},
               {Spin::beta /* 9*/,Ethane::nuclei.atoms().positionsVector()[7]},
               // C-C Bond
               {Spin::alpha/*10*/,inbetween(Ethane::nuclei.atoms(),{0,1},0.25)},
               {Spin::beta /*11*/,inbetween(Ethane::nuclei.atoms(),{1,0},0.25)},
               // C-H Bonds
               {Spin::beta /*12*/,inbetween(Ethane::nuclei.atoms(),{0,2},0.25)},
               {Spin::beta /*13*/,Ethane::nuclei.atoms().positionsVector()[3]},
               {Spin::beta /*14*/,inbetween(Ethane::nuclei.atoms(),{0,4},0.25)},
               {Spin::alpha/*15*/,inbetween(Ethane::nuclei.atoms(),{1,5},0.25)},
               {Spin::alpha/*16*/,Ethane::nuclei.atoms().positionsVector()[6]},
               {Spin::alpha/*17*/,inbetween(Ethane::nuclei.atoms(),{1,7},0.25)}
            })};

    std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permIndices[0] << 0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15; // 120° rotation around z
    permIndices[1] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}