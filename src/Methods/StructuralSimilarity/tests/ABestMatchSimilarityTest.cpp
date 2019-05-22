//
// Created by heuer on 03.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <limits>
#include <random>
#include "Hungarian.h"

#include <iomanip>

using namespace testing;
using namespace SOAP;

class ABestMatchSimilarityTest : public ::testing::Test {
public:
    MolecularGeometry A, B, C, D;

    double simRadius, soapThreshold;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);

        simRadius = 0.1;
        soapThreshold = 1.0;
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

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, simRadius, soapThreshold);

    Eigen::VectorXi a(6),b(6),c(6),d(6),e(6),f(6),g(6),h(6);

    a << 0,1,2,3,4,5;
    b << 0,1,2,5,4,3;
    c << 0,1,3,5,4,2;
    d << 2,1,0,3,4,5;
    e << 2,1,3,5,4,0;
    f << 5,1,2,0,4,3;
    g << 5,1,3,0,4,2;
    h << 5,1,3,2,4,0;

    for (auto& i : results) {
        ASSERT_THAT(i.permutation.indices(), AnyOf(a,b,c,d,e,f,g,h));
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

    auto results = BestMatch::SOAPSimilarity::getAllBestMatchResults(specA,specB, simRadius, soapThreshold);

    Eigen::VectorXi a(6),b(6),c(6),d(6),e(6),f(6),g(6),h(6);

    a << 0,1,2,3,4,5;
    b << 0,1,2,5,4,3;
    c << 0,1,3,5,4,2;
    d << 2,1,0,3,4,5;
    e << 2,1,3,5,4,0;
    f << 5,1,2,0,4,3;
    g << 5,1,3,0,4,2;
    h << 5,1,3,2,4,0;

    for (auto& i : results) {
        ASSERT_THAT(i.permutation.indices(), AnyOf(a,b,c,d,e,f,g,h));
    }
}
