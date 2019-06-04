//
// Created by heuer on 03.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <limits>
#include <random>
#include <Metrics.h>
#include <Hungarian.h>

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


TEST_F(ABestMatchSimilarityTest, VarySimilarEnvironments) {

}

/*
TEST_F(ABestMatchSimilarityTest, IndexReordering) {
    std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> test =
            {
            {{{{0,0},{1,1},{4,4},{2,2}}}},
            {{{7,7},{8,8}}},
            {{{5,5},{3,3}}}
            };

    auto result = BestMatch::SOAPSimilarity::obtainIndexReorderingPermutationOverAllBlocks(test);

    ASSERT_THAT(result, ElementsAre(0,1,2,6,3,7,4,5));
}*/

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    ParticleKit::create(A);


    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_Shaked) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    ParticleKit::create(A);

    while(Metrics::positionalNormsVectorNorm<Eigen::Dynamic,2>(
            A.electrons().positionsVector(),
            B.electrons().positionsVector()) == 0.0)
        A.electrons().positionsVector().shake(distanceTolerance/2.0);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, shakeSoapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_LT(i.metric, 1.0);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Alchemical) {
    General::settings.mode = General::Mode::alchemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    ParticleKit::create(A);


    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Alchemical_Shaked) {
    General::settings.mode = General::Mode::alchemical;

    auto B = TestMolecules::H4::fourAlpha;
    auto A = B;
    ParticleKit::create(A);

    while(Metrics::positionalNormsVectorNorm<Eigen::Dynamic,2>(
            A.electrons().positionsVector(),
            B.electrons().positionsVector()) == 0.0)
        A.electrons().positionsVector().shake(distanceTolerance/2.0);

    // Permute electrons in A
    //Eigen::PermutationMatrix<Eigen::Dynamic> perm(A.electrons().numberOfEntities());
    //perm.setIdentity();
    ////auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));
    //auto rng = std::default_random_engine(static_cast<unsigned long>(0));
    //std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
    //A.electrons().permute(perm);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, shakeSoapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(4,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_LT(i.metric, 1.0);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_BoraneIonic) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::BH3::ionicMinimal;
    auto A = TestMolecules::BH3::ionicMinimalRotated;
    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3;
    permsIndices[1] << 1,0, 3,2;

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, soapThreshold);

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        std::cout << i.metric << " ,  " << i.permutation.indices().transpose() << std::endl;
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_EthaneIonic) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::Ethane::doublyIonicMinimal;
    auto A = TestMolecules::Ethane::doublyIonicMinimalDoublyRotated;
    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3, 4,5, 6,7;
    permsIndices[1] << 1,0, 3,2, 5,4, 7,6;
    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, soapThreshold);

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

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

    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    std::vector<Eigen::VectorXi> permsIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1, 2,3;
    permsIndices[1] << 2,3, 0,1;
    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(specA, specB, distanceTolerance, soapThreshold);

    std::wcerr << "check. why the other correct permutations are not found " << std::endl;
    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        std::cout << i.metric << " ,  " << i.permutation.indices().transpose() << std::endl;
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
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

    ParticleKit::create(A);
    auto permutee = MolecularSpectrum(A);
    auto reference = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(
            permutee, reference, distanceTolerance, soapThreshold);


    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    expectedPermIndices[0] << 2,0,1,5,3,4; // 120° rotation around z
    expectedPermIndices[1] << 0,2,1,3,5,4; // reflection along H2-C0-C1-H5 plane

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(expectedPermIndices));
    }
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

    ParticleKit::create(A);
    auto permutee = MolecularSpectrum(A);
    auto reference = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,2,3,  6,4,5, 9,7,8,  10,11,  14,12,13, 17,15,16; // 120° rotation around z
    expectedPermIndices[1] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(expectedPermIndices));
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

    ParticleKit::create(A);
    auto permutee = MolecularSpectrum(A);
    auto reference = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(
            permutee, reference, distanceTolerance, soapThreshold);


    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    expectedPermIndices[0] << 2,0,1,5,3,4; // 120° rotation around z
    expectedPermIndices[1] << 0,2,1,3,5,4; // reflection along H2-C0-C1-H5 plane

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(expectedPermIndices));
    }
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

    ParticleKit::create(A);
    auto permutee = MolecularSpectrum(A);
    auto reference = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,2,3,  6,4,5, 9,7,8,  10,11,  14,12,13, 17,15,16; // 120° rotation around z
    expectedPermIndices[1] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane

    ASSERT_EQ(results.size(), expectedPermIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(expectedPermIndices));
    }
}