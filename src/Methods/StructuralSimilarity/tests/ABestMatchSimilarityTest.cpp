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

using namespace testing;
using namespace SOAP;

class ABestMatchSimilarityTest : public ::testing::Test {
public:
    MolecularGeometry A, B, C, D;

    double distanceTolerance, soapThreshold, shakeSoapThreshold, eps;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);

        distanceTolerance = 0.1;
        soapThreshold = 1.0;
        shakeSoapThreshold = 0.90;
        eps = std::numeric_limits<double>::epsilon()*10;
    };
};


TEST_F(ABestMatchSimilarityTest, PrermuteEnvironmentsToLabSystem) {

    auto b = TestMolecules::BH3::ionicMinimal;
    auto a = TestMolecules::BH3::ionicMinimalRotatedPermuted;

    ParticleKit::create(a);
    General::settings.mode = General::Mode::chemical;

    Radial::settings.nmax = 2;
    Angular::settings.lmax = 2;

    auto specA = MolecularSpectrum(a);
    auto specB = MolecularSpectrum(b);

    auto environmentalSimilarities = BestMatch::SOAPSimilarity::calculateEnvironmentalSimilarityMatrix(specA,specB);

    auto fromKitA = ParticleKit::fromKitPermutation(a.electrons());
    auto fromKitB = ParticleKit::fromKitPermutation(b.electrons());

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

    auto result = BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(mat,1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_THAT(result[0], ElementsAre(0,2));
    ASSERT_THAT(result[1], ElementsAre(1,3));
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndices2) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,1,0,0,\
    1,1,0,0,\
    0,0,1,1,\
    0,0,1,1;

    auto result = BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(mat,1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_THAT(result[0], ElementsAre(0,1));
    ASSERT_THAT(result[1], ElementsAre(2,3));
}

TEST_F(ABestMatchSimilarityTest, IndexReordering) {
    std::deque<std::vector<std::deque<Eigen::Index>>> test = {{{0,1,4,2}},{{7,8}},{{5,3}}};

    auto result = BestMatch::SOAPSimilarity::obtainIndexReorderingPermutationOverAllBlocks(test);

    ASSERT_THAT(result, ElementsAre(0,1,2,6,3,7,4,5));
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical) {
    auto B = TestMolecules::H4::sixElectrons;
    auto A = B;

    ParticleKit::create(A);
    General::settings.mode = General::Mode::chemical;

    Radial::settings.nmax = 2;
    Angular::settings.lmax = 2;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(8,Eigen::VectorXi(6));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,1,2,5,4,3;
    permsIndices[2] << 0,1,3,5,4,2;
    permsIndices[3] << 2,1,0,3,4,5;
    permsIndices[4] << 2,1,3,5,4,0;
    permsIndices[5] << 5,1,2,0,4,3;
    permsIndices[6] << 5,1,3,0,4,2;
    permsIndices[7] << 5,1,3,2,4,0;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Chemical_Shaked) {
    auto B = TestMolecules::H4::sixElectrons;
    auto A = B;

    while(Metrics::positionalNormsVectorNorm<Eigen::Dynamic,2>(
            A.electrons().positionsVector(),
            B.electrons().positionsVector()) == 0.0)
        A.electrons().positionsVector().shake(distanceTolerance/2.0);

    ParticleKit::create(A);
    General::settings.mode = General::Mode::chemical;

    Radial::settings.nmax = 3;
    Radial::settings.sigmaAtom = 2.0;
    Angular::settings.lmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 1.0;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, distanceTolerance, shakeSoapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(8,Eigen::VectorXi(6));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,1,2,5,4,3;
    permsIndices[2] << 0,1,3,5,4,2;
    permsIndices[3] << 2,1,0,3,4,5;
    permsIndices[4] << 2,1,3,5,4,0;
    permsIndices[5] << 5,1,2,0,4,3;
    permsIndices[6] << 5,1,3,0,4,2;
    permsIndices[7] << 5,1,3,2,4,0;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_LT(i.metric, 1.0);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Alchemical) {
    auto B = TestMolecules::H4::sixElectrons;
    auto A = B;

    ParticleKit::create(A);
    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;


    Radial::settings.nmax = 2;
    Angular::settings.lmax = 2;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, distanceTolerance, soapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(8,Eigen::VectorXi(6));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,1,2,5,4,3;
    permsIndices[2] << 0,1,3,5,4,2;
    permsIndices[3] << 2,1,0,3,4,5;
    permsIndices[4] << 2,1,3,5,4,0;
    permsIndices[5] << 5,1,2,0,4,3;
    permsIndices[6] << 5,1,3,0,4,2;
    permsIndices[7] << 5,1,3,2,4,0;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_NEAR(i.metric, 1.0, eps);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}

TEST_F(ABestMatchSimilarityTest, FindDistanceConservingPermutations_Alchemical_Shaked) {
    auto B = TestMolecules::H4::sixElectrons;
    auto A = B;

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

    ParticleKit::create(A);
    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;

    Radial::settings.nmax = 3;
    Radial::settings.sigmaAtom = 2.0;
    Angular::settings.lmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 1.0;
    //Cutoff::settings.centerWeight = 0.1;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, distanceTolerance, shakeSoapThreshold);

    std::vector<Eigen::VectorXi> permsIndices(8,Eigen::VectorXi(6));
    permsIndices[0] << 0,1,2,3,4,5;
    permsIndices[1] << 0,1,2,5,4,3;
    permsIndices[2] << 0,1,3,5,4,2;
    permsIndices[3] << 2,1,0,3,4,5;
    permsIndices[4] << 2,1,3,5,4,0;
    permsIndices[5] << 5,1,2,0,4,3;
    permsIndices[6] << 5,1,3,0,4,2;
    permsIndices[7] << 5,1,3,2,4,0;

    ASSERT_EQ(results.size(), permsIndices.size());
    for (auto& i : results) {
        ASSERT_LT(i.metric, 1.0);
        ASSERT_THAT(i.permutation.indices(), AnyOfArray(permsIndices));
    }
}