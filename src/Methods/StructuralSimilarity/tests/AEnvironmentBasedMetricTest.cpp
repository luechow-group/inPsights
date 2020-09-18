// Copyright (C) 2019-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <limits>
#include <random>
#include <EnvironmentBasedMetric.h>
#include <StructuralSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <Metrics.h>
#include <Combinatorics.h>
#include <Enumerate.h>

using namespace testing;
using namespace SOAP;

class AEnvironmentBasedMetricTest : public ::testing::Test {
public:
    using Pair = std::pair<Eigen::Index,Eigen::Index>;
    using MolPermList = std::vector<std::pair<std::vector<int>,std::vector<int>>>;
    double distanceTolerance, soapThreshold, shakeSoapThreshold;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        //spdlog::set_level(spdlog::level::debug);

        distanceTolerance = 0.1;
        soapThreshold = 1.0;
        shakeSoapThreshold = 0.90;

        Radial::settings.nmax = 2;
        Radial::settings.sigmaAtom = 2.0;
        Angular::settings.lmax = 2;
        Cutoff::settings.radius = 6.0;
        Cutoff::settings.width = 0.0;
        General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    };

    void routine(const MolecularGeometry &A, const MolecularGeometry &B,
                 const MolPermList &expectedPermutationIndices,
                 double distTolerance, double soapThresh, bool greaterThan = false
                 ) {

        double eps = SOAP::General::settings.comparisonEpsilon();

        ParticleKit::create(A);
        ASSERT_TRUE(ParticleKit::isSubsetQ(A));
        ASSERT_TRUE(ParticleKit::isSubsetQ(B));

        auto specA = MolecularSpectrum(A);
        auto specB = MolecularSpectrum(B);

        auto results = Metrics::Similarity::SOAPBased::getBestMatchResults(
                specA, specB,
                distTolerance, soapThresh, eps);

        std::sort(results.begin(), results.end());

        spdlog::debug("Expected:");
        for (auto[i, expected] : enumerate(expectedPermutationIndices))
            spdlog::debug("{}: metric {} {}, nuclear = {}, electronic = {} (lab system)",
                    i, greaterThan ? ">" : "=", soapThresh,
                    ToString::stdvectorIntToString(expected.first),
                    ToString::stdvectorIntToString(expected.second));


        unsigned removedCounter = 0;
        for(auto it = results.begin(); it != results.end(); it++) {
            if(it->metric < (soapThresh-eps)) {
                results.erase(it--);
                removedCounter++;
            }
        }
        spdlog::debug("Removed {} result(s) below threshold.", removedCounter);

        spdlog::debug("Got:");
        for (auto[i, result] : enumerate(results)) {
            auto molecularPermResult = A.splitAllParticlePermutation(result.permutation);

            spdlog::debug("{}: metric = {:01.16f}, nuclear = {}, electronic = {} (lab system)", i, result.metric,
                          ToString::vectorXiToString(molecularPermResult.nuclearPermutation.indices()),
                          ToString::vectorXiToString(molecularPermResult.electronicPermutation.indices()));
        }

        size_t i = 0;
        while (i < results.size() || i < expectedPermutationIndices.size()) {

            if (i >= results.size()) {
                spdlog::debug(" {}: missing", i);
                i++;
                continue;
            }

            if (i < expectedPermutationIndices.size()) {
                if (!greaterThan) ASSERT_NEAR(results[i].metric, soapThresh, eps);
                else ASSERT_GT(results[i].metric, soapThresh);

                auto molPerm = A.splitAllParticlePermutation(results[i].permutation);
                auto nucIndices = molPerm.nuclearPermutation.indices();
                std::vector<int> nucIndicesAsVector(nucIndices.data(), nucIndices.data()+ nucIndices.size());
                ASSERT_EQ(nucIndicesAsVector, expectedPermutationIndices[i].first);

                auto elecIndices = molPerm.electronicPermutation.indices();
                std::vector<int> elecIndicesAsVector(elecIndices.data(), elecIndices.data()+ elecIndices.size());
                ASSERT_EQ(elecIndicesAsVector,expectedPermutationIndices[i].second);
            } else {
                spdlog::debug("Found unexpected result {}", i);
            };
            i++;

        }
        ASSERT_EQ(results.size(), expectedPermutationIndices.size());
    }
};

TEST_F(AEnvironmentBasedMetricTest, PrermuteEnvironmentsToLabSystem) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::BH3::ionicMinimal;
    auto A = TestMolecules::BH3::ionicMinimalRotatedPermuted;
    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto environmentalSimilarities = SOAP::StructuralSimilarity::correlationMatrix(specA, specB);

    auto fromKitA = ParticleKit::fromKitPermutation(A);
    auto fromKitB = ParticleKit::fromKitPermutation(B);

    auto environmentalSimilaritiesInLabSystem = fromKitA*environmentalSimilarities*fromKitB;

    auto eps = std::numeric_limits<double>::epsilon()*10;

    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(0,2), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(1,0), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(1,1), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(2,0), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(2,1), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(3,3), 1.0, eps);

    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(4,4), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(4,5), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(5,6), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(5,7), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(6,4), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(6,5), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(7,6), 1.0, eps);
    ASSERT_NEAR(environmentalSimilaritiesInLabSystem(7,7), 1.0, eps);
}

TEST_F(AEnvironmentBasedMetricTest, GrowingPerm) {

    auto growingPerm = Metrics::Similarity::SOAPBased::GrowingPerm({0, 1, 2}, {});

    ASSERT_TRUE(growingPerm.add({0,1}));
    ASSERT_THAT(growingPerm.remainingPermuteeIndices_,ElementsAre(1,2));

    ASSERT_EQ(growingPerm.chainOfSwaps_[0].first, 0);
    ASSERT_EQ(growingPerm.chainOfSwaps_[0].second, 1);

    ASSERT_FALSE(growingPerm.add({0,2}));
}

TEST_F(AEnvironmentBasedMetricTest, FindPossiblePermutations1) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,0,1,0,\
    0,1,0,1,\
    1,0,1,0,\
    0,1,0,1;

    auto matches = Metrics::Similarity::SOAPBased::findEnvironmentMatches(mat, 1.0 );
    ASSERT_THAT(matches[0].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(matches[1].permuteeIndices, ElementsAre(1, 3));
    ASSERT_THAT(matches[2].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(matches[3].permuteeIndices, ElementsAre(1, 3));
    for (size_t i = 0; i < matches.size(); ++i)
        ASSERT_EQ(matches[i].referenceIndex, i);

    auto dependentMatches = Metrics::Similarity::SOAPBased::clusterDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(dependentMatches[0][1].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(dependentMatches[1][0].permuteeIndices, ElementsAre(1, 3));
    ASSERT_THAT(dependentMatches[1][1].permuteeIndices, ElementsAre(1, 3));

    auto possiblePermutations0 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[0]);
    ASSERT_EQ(possiblePermutations0.size(), 2);
    ASSERT_THAT(possiblePermutations0[0].chainOfSwaps_, ElementsAre(Pair(0,0),Pair(2,2)));
    ASSERT_THAT(possiblePermutations0[1].chainOfSwaps_, ElementsAre(Pair(2,0),Pair(0,2)));
    for (auto & i : possiblePermutations0)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());

    auto possiblePermutations1 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[1]);
    ASSERT_EQ(possiblePermutations1.size(), 2);
    ASSERT_THAT(possiblePermutations1[0].chainOfSwaps_, ElementsAre(Pair(1,1),Pair(3,3)));
    ASSERT_THAT(possiblePermutations1[1].chainOfSwaps_, ElementsAre(Pair(3,1),Pair(1,3)));
    for (auto & i : possiblePermutations1)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());
}

TEST_F(AEnvironmentBasedMetricTest, FindPossiblePermutations2) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,1,0,1,\
    1,1,0,0,\
    0,0,1,0,\
    1,0,0,1;

    auto matches = Metrics::Similarity::SOAPBased::findEnvironmentMatches(mat, 1.0);
    ASSERT_THAT(matches[0].permuteeIndices, ElementsAre(0, 1, 3));
    ASSERT_THAT(matches[1].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(matches[2].permuteeIndices, ElementsAre(2));
    ASSERT_THAT(matches[3].permuteeIndices, ElementsAre(0, 3));
    for (size_t i = 0; i < matches.size(); ++i)
        ASSERT_EQ(matches[i].referenceIndex, i);

    auto dependentMatches = Metrics::Similarity::SOAPBased::clusterDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeIndices, ElementsAre(0, 1, 3));
    ASSERT_THAT(dependentMatches[0][1].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(dependentMatches[0][2].permuteeIndices, ElementsAre(0, 3));

    ASSERT_THAT(dependentMatches[1][0].permuteeIndices, ElementsAre(2));
    ASSERT_EQ(dependentMatches[1][0].referenceIndex, 2);

    auto possiblePermutations0 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[0]);
    ASSERT_EQ(possiblePermutations0.size(), 3);
    ASSERT_THAT(possiblePermutations0[0].chainOfSwaps_, ElementsAre(Pair(0,0),Pair(1,1),Pair(3,3)));
    ASSERT_THAT(possiblePermutations0[1].chainOfSwaps_, ElementsAre(Pair(1,0),Pair(0,1),Pair(3,3)));
    ASSERT_THAT(possiblePermutations0[2].chainOfSwaps_, ElementsAre(Pair(3,0),Pair(1,1),Pair(0,3)));
    for (auto & i : possiblePermutations0)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());

    auto possiblePermutations1 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[1]);
    ASSERT_EQ(possiblePermutations1.size(), 1);
    ASSERT_THAT(possiblePermutations1[0].chainOfSwaps_, ElementsAre(Pair(2,2)));
    for (auto & i : possiblePermutations1)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());
}

TEST_F(AEnvironmentBasedMetricTest, FindPossiblePermutations3) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,1,0,0,\
    1,1,0,0,\
    0,0,1,1,\
    0,0,1,1;

    auto matches = Metrics::Similarity::SOAPBased::findEnvironmentMatches(mat, 1.0);
    ASSERT_THAT(matches[0].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(matches[1].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(matches[2].permuteeIndices, ElementsAre(2, 3));
    ASSERT_THAT(matches[3].permuteeIndices, ElementsAre(2, 3));
    for (size_t i = 0; i < matches.size(); ++i)
        ASSERT_EQ(matches[i].referenceIndex, i);

    auto dependentMatches = Metrics::Similarity::SOAPBased::clusterDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(dependentMatches[0][1].permuteeIndices, ElementsAre(0, 1));
    ASSERT_THAT(dependentMatches[1][0].permuteeIndices, ElementsAre(2, 3));
    ASSERT_THAT(dependentMatches[1][1].permuteeIndices, ElementsAre(2, 3));

    auto possiblePermutations0 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[0]);
    ASSERT_EQ(possiblePermutations0.size(), 2);
    ASSERT_THAT(possiblePermutations0[0].chainOfSwaps_, ElementsAre(Pair(0,0),Pair(1,1)));
    ASSERT_THAT(possiblePermutations0[1].chainOfSwaps_, ElementsAre(Pair(1,0),Pair(0,1)));
    for (const auto & i : possiblePermutations0)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());

    auto possiblePermutations1 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[1]);
    ASSERT_EQ(possiblePermutations1.size(), 2);
    ASSERT_THAT(possiblePermutations1[0].chainOfSwaps_, ElementsAre(Pair(2,2),Pair(3,3)));
    ASSERT_THAT(possiblePermutations1[1].chainOfSwaps_, ElementsAre(Pair(3,2),Pair(2,3)));
    for (const auto & i : possiblePermutations1)
        ASSERT_TRUE(i.remainingPermuteeIndices_.empty());
}

TEST_F(AEnvironmentBasedMetricTest, FindPossiblePermutations4) {
    Eigen::MatrixXd mat(6,6);
    mat <<
    1,1,0,0,0,0,\
    0,0,1,0,0,0,\
    1,1,0,0,0,0,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1;

    auto matches = Metrics::Similarity::SOAPBased::findEnvironmentMatches(mat, 1.0);
    ASSERT_THAT(matches[0].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(matches[1].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(matches[2].permuteeIndices, ElementsAre(1));
    ASSERT_THAT(matches[3].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(matches[4].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(matches[5].permuteeIndices, ElementsAre(3, 4, 5));
    for (size_t i = 0; i < matches.size(); ++i)
        ASSERT_EQ(matches[i].referenceIndex, i);

    auto dependentMatches = Metrics::Similarity::SOAPBased::clusterDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeIndices, ElementsAre(0, 2));
    ASSERT_THAT(dependentMatches[0][1].permuteeIndices, ElementsAre(0, 2));

    ASSERT_THAT(dependentMatches[1][0].permuteeIndices, ElementsAre(1));

    ASSERT_THAT(dependentMatches[2][0].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(dependentMatches[2][1].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(dependentMatches[2][2].permuteeIndices, ElementsAre(3, 4, 5));

    auto possiblePermutations0 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[0]);
    ASSERT_EQ(possiblePermutations0.size(), 2);
    ASSERT_THAT(possiblePermutations0[0].chainOfSwaps_, ElementsAre(Pair(0,0),Pair(2,1)));
    ASSERT_THAT(possiblePermutations0[1].chainOfSwaps_, ElementsAre(Pair(2,0),Pair(0,1)));

    auto possiblePermutations1 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[1]);
    ASSERT_EQ(possiblePermutations1.size(), 1);
    ASSERT_THAT(possiblePermutations1[0].chainOfSwaps_, ElementsAre(Pair(1,2)));

    auto possiblePermutations2 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[2]);
    ASSERT_EQ(possiblePermutations2.size(), 6);
    ASSERT_THAT(possiblePermutations2[0].chainOfSwaps_, ElementsAre(Pair(3,3),Pair(4,4),Pair(5,5)));
    ASSERT_THAT(possiblePermutations2[1].chainOfSwaps_, ElementsAre(Pair(3,3),Pair(5,4),Pair(4,5)));
    ASSERT_THAT(possiblePermutations2[2].chainOfSwaps_, ElementsAre(Pair(4,3),Pair(3,4),Pair(5,5)));
    ASSERT_THAT(possiblePermutations2[3].chainOfSwaps_, ElementsAre(Pair(4,3),Pair(5,4),Pair(3,5)));
    ASSERT_THAT(possiblePermutations2[4].chainOfSwaps_, ElementsAre(Pair(5,3),Pair(3,4),Pair(4,5)));
    ASSERT_THAT(possiblePermutations2[5].chainOfSwaps_, ElementsAre(Pair(5,3),Pair(4,4),Pair(3,5)));
}

TEST_F(AEnvironmentBasedMetricTest, FindPossiblePermutations5) {
    Eigen::MatrixXd mat(6,6);
    mat <<
    0,1,0,0,0,0,\
    0,0,1,0,0,0,\
    1,0,0,0,0,0,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1;

    auto matches = Metrics::Similarity::SOAPBased::findEnvironmentMatches(mat, 1.0);
    ASSERT_THAT(matches[0].permuteeIndices, ElementsAre(2));
    ASSERT_THAT(matches[1].permuteeIndices, ElementsAre(0));
    ASSERT_THAT(matches[2].permuteeIndices, ElementsAre(1));
    ASSERT_THAT(matches[3].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(matches[4].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(matches[5].permuteeIndices, ElementsAre(3, 4, 5));
    for (size_t i = 0; i < matches.size(); ++i)
        ASSERT_EQ(matches[i].referenceIndex, i);

    auto dependentMatches = Metrics::Similarity::SOAPBased::clusterDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeIndices, ElementsAre(2));

    ASSERT_THAT(dependentMatches[1][0].permuteeIndices, ElementsAre(0));

    ASSERT_THAT(dependentMatches[2][0].permuteeIndices, ElementsAre(1));

    ASSERT_THAT(dependentMatches[3][0].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(dependentMatches[3][1].permuteeIndices, ElementsAre(3, 4, 5));
    ASSERT_THAT(dependentMatches[3][2].permuteeIndices, ElementsAre(3, 4, 5));

    auto possiblePermutations0 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[0]);
    ASSERT_EQ(possiblePermutations0.size(), 1);
    ASSERT_THAT(possiblePermutations0[0].chainOfSwaps_, ElementsAre(Pair(2,0)));

    auto possiblePermutations1 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[1]);
    ASSERT_EQ(possiblePermutations1.size(), 1);
    ASSERT_THAT(possiblePermutations1[0].chainOfSwaps_, ElementsAre(Pair(0,1)));

    auto possiblePermutations2 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[2]);
    ASSERT_EQ(possiblePermutations2.size(), 1);
    ASSERT_THAT(possiblePermutations2[0].chainOfSwaps_, ElementsAre(Pair(1,2)));

    auto possiblePermutations3 = Metrics::Similarity::SOAPBased::findPossiblePermutations(dependentMatches[3]);
    ASSERT_EQ(possiblePermutations3.size(), 6);
    ASSERT_THAT(possiblePermutations3[0].chainOfSwaps_, ElementsAre(Pair(3,3),Pair(4,4),Pair(5,5)));
    ASSERT_THAT(possiblePermutations3[1].chainOfSwaps_, ElementsAre(Pair(3,3),Pair(5,4),Pair(4,5)));
    ASSERT_THAT(possiblePermutations3[2].chainOfSwaps_, ElementsAre(Pair(4,3),Pair(3,4),Pair(5,5)));
    ASSERT_THAT(possiblePermutations3[3].chainOfSwaps_, ElementsAre(Pair(4,3),Pair(5,4),Pair(3,5)));
    ASSERT_THAT(possiblePermutations3[4].chainOfSwaps_, ElementsAre(Pair(5,3),Pair(3,4),Pair(4,5)));
    ASSERT_THAT(possiblePermutations3[5].chainOfSwaps_, ElementsAre(Pair(5,3),Pair(4,4),Pair(3,5)));
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_alchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;

    MolPermList expectedPerms{
            {{3,2,1,0},{1,2,3,0}},
            {{3,2,1,0},{3,2,1,0}}
    };
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_ionic_reflected) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflected;

    MolPermList expectedPerms{
            {{3,2,1,0},{0,1,2,3}},
            {{3,2,1,0},{2,1,0,3}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPerms = {expectedPerms[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_ionic_reflected_alpha_permuted) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedAlphaPermuted;

    MolPermList expectedPerms{
            {{3,2,1,0},{1,0,2,3}},
            {{3,2,1,0},{2,0,1,3}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPerms = {expectedPerms[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_ionic_reflected_beta_permuted) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedBetaPermuted;

    MolPermList expectedPerms{
            {{3,2,1,0},{0,1,3,2}},
            {{3,2,1,0},{3,1,0,2}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPerms = {expectedPerms[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_ionic_reflected_reordered_numbering) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering;

    MolPermList expectedPerms{
            {{3,2,1,0},{1,2,3,0}},
            {{3,2,1,0},{3,2,1,0}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    // both perms are not conserving in chemical mode
    expectedPerms = {};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_not_in_particle_kit) {
    auto B = TestMolecules::H4::linear::ionicA;
    auto A = TestMolecules::H4::linear::ionicNotinParticleKitSystem;

    MolPermList expectedPerms{
            {{0,1,2,3},{0,1,2,3}},
            {{0,1,2,3},{2,1,0,3}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPerms = {expectedPerms[1]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_real_maxima1) {
    auto A = TestMolecules::H4::linear::ionicRealMax1Var1;
    auto B = TestMolecules::H4::linear::ionicRealMax1Var2;

    MolPermList expectedPerms{
            {{3,2,1,0},{0,2,1,3}},
            {{3,2,1,0},{3,2,1,0}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    General::settings.mode = General::Mode::chemical;
    auto envSimMat = SOAP::StructuralSimilarity::correlationMatrix(MolecularSpectrum(A), MolecularSpectrum(B));
    ASSERT_GT(envSimMat.minCoeff(), 0.0);
    ASSERT_LT(envSimMat.maxCoeff(), 1.0); // no equivalent environments exist in chemical mode.
    expectedPerms = {};
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4linear_real_maxima2) {
    auto A = TestMolecules::H4::linear::ionicRealMax2Var1;
    auto B = TestMolecules::H4::linear::ionicRealMax2Var2;

    MolPermList expectedPerms{
            {{3,2,1,0},{1,0,3,2}},
            {{3,2,1,0},{1,3,0,2}},
            {{3,2,1,0},{2,0,3,1}},
            {{3,2,1,0},{2,3,0,1}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);

    expectedPerms = {expectedPerms[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPerms, distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4ring) {
    auto B = TestMolecules::H4::ring::fourAlpha;
    auto A = B;

    MolPermList expectedPerms{
            {{0,1,2,3},{0,1,2,3}},
            {{0,1,3,2},{0,1,3,2}},
            {{1,0,2,3},{1,0,2,3}},
            {{1,0,3,2},{1,0,3,2}},
            {{2,3,0,1},{2,3,0,1}},
            {{2,3,1,0},{2,3,1,0}},
            {{3,2,0,1},{3,2,0,1}},
            {{3,2,1,0},{3,2,1,0}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);

    General::settings.mode = General::Mode::chemical;
    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, H4ring_Shaked) {
    auto B = TestMolecules::H4::ring::fourAlpha;
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

        MolPermList expectedPerms{
                {{0,1,2,3},{0,1,2,3}},
                {{0,1,3,2},{0,1,3,2}},
                {{1,0,2,3},{1,0,2,3}},
                {{1,0,3,2},{1,0,3,2}},
                {{2,3,0,1},{2,3,0,1}},
                {{2,3,1,0},{2,3,1,0}},
                {{3,2,0,1},{3,2,0,1}},
                {{3,2,1,0},{3,2,1,0}}
        };

        General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
        General::settings.mode = General::Mode::alchemical;
        routine(A, B, expectedPerms, distanceTolerance, shakeSoapThreshold, true);

        General::settings.mode = General::Mode::chemical;
        routine(A, B, expectedPerms, distanceTolerance, shakeSoapThreshold, true);
    }
}

TEST_F(AEnvironmentBasedMetricTest, FindDistanceConservingPermutations_Chemical_BoraneLEO1) {
    General::settings.mode = General::Mode::chemical;

    auto nuclei = TestMolecules::BH3::nuclei;
    const MolecularGeometry A = {
            TestMolecules::BH3::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                {Spin::alpha, nuclei.atoms().positionsVector()[3]}
            })};

    MolPermList expectedPerms{
            {{0,1,2,3},{0,1,2}},
            {{0,1,3,2},{0,2,1}},
            {{0,2,1,3},{1,0,2}},
            {{0,2,3,1},{1,2,0}},
            {{0,3,1,2},{2,0,1}},
            {{0,3,2,1},{2,1,0}}
    };

    routine(A,A,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, BH3Covalent_Chemical) {
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

    MolPermList expectedPerms{
            {{0,1,2,3},{0,1,2,3,4,5}},
            {{0,1,3,2},{0,2,1,3,5,4}},
            {{0,2,1,3},{1,0,2,4,3,5}},
            {{0,2,3,1},{1,2,0,4,5,3}},
            {{0,3,1,2},{2,0,1,5,3,4}},
            {{0,3,2,1},{2,1,0,5,4,3}}
    };

    routine(A,A,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, BH3CovalentPermuted_Chemical) {
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

    MolPermList expectedPerms{
            {{0,1,2,3},{0,1,2,3,4,5}},
            {{0,1,3,2},{0,4,2,5,1,3}},
            {{0,2,1,3},{1,0,3,2,4,5}},
            {{0,2,3,1},{1,4,3,5,0,2}},
            {{0,3,1,2},{4,0,5,2,1,3}},
            {{0,3,2,1},{4,1,5,3,0,2}}
    };

    routine(A,A,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, BH3_ThreeIndependentUnsimilarEnvironments_Chemical) {
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

    MolPermList expectedPerms{
            {{0,2,1,3},{1,0,2}}
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}


TEST_F(AEnvironmentBasedMetricTest, BH3Ionic_Chemical) {
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

    MolPermList expectedPerms{
            {{0, 2,3, 1},{0,1, 2,3}},
            {{0, 3,2, 1},{1,0, 3,2}}
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, EthaneDoublyIonicMinimal) {
    General::settings.mode = General::Mode::chemical;

    MolecularGeometry A = {
            TestMolecules::Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[6]},
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[3]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[3]}}
            )};

    MolecularGeometry B = {
            TestMolecules::Ethane::nuclei.atoms(),
            ElectronsVector({
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[2]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[2]},
                {Spin::alpha,TestMolecules::Ethane::nuclei.atoms().positionsVector()[5]},
                {Spin::beta, TestMolecules::Ethane::nuclei.atoms().positionsVector()[5]}}
            )};


    MolPermList expectedPerms{
            {{0,1, 3,2,4, 6,5,7}, {2,3, 0,1}}, // reflection along H4-C0-C1-C7 Plane
            {{0,1, 4,2,3, 7,5,6}, {2,3, 0,1}}, // 120° rotation around z
            {{1,0, 6,5,7, 3,2,4}, {0,1, 2,3}}, // reflection in xy + 60° rotation around z
            {{1,0, 7,5,6, 4,2,3}, {0,1, 2,3}}, // reflection in xy - 60° rotation around z
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, EthaneSinglyIonicMinimal) {
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

    MolPermList expectedPerms{
            {{0,1, 2,4,3, 5,7,6},{0, 2, 1, 3, 5, 4}},// reflection along H2-C0-C1-H5 plane
            {{0,1, 3,4,2, 6,7,5},{1, 2, 0, 4, 5, 3}} // 120° rotation around z
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, EthaneSinglyIonicMinimal_shaked_alchemical) {
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

    MolPermList expectedPerms{
            {{0,1, 2,4,3, 5,7,6},{0, 2, 1, 3, 5, 4}},// reflection along H2-C0-C1-H5 plane
            {{0,1, 3,4,2, 6,7,5},{1, 2, 0, 4, 5, 3}} // 120° rotation around z
    };

    std::random_device rd;
    std::uniform_int_distribution<unsigned long> dist(0, 9999);
    auto randomSeed = dist(rd);
    std::cout << "random seed: " << randomSeed << std::endl;

    auto rng = std::default_random_engine(randomSeed);
    A.electrons().positionsVector().shake(distanceTolerance / 10.0, rng);

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    SOAP::General::settings.comparisonEpsilon = 1e-8;

    routine(A, B, expectedPerms, distanceTolerance, shakeSoapThreshold, true);
}

TEST_F(AEnvironmentBasedMetricTest, EthaneSinglyIonic) {
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

    MolPermList expectedPerms{
            {{0,1, 2,4,3, 5,7,6},{0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16}},
            {{0,1, 3,4,2, 6,7,5},{0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15}}
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}


TEST_F(AEnvironmentBasedMetricTest, EthaneSinglyIonic10RandomIndexSwapPermutations) {

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

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;

    std::vector<int> indices(A.electrons().numberOfEntities());
    std::iota(indices.begin(), indices.end(), 0);

    auto allIndexSwaps = Combinatorics::Permutations<int>(indices).get();

    std::random_device random_device;
    auto randomSeed = random_device();
    std::cout << "Random seed: " << randomSeed << std::endl;
    std::mt19937 engine{random_device()};
    std::uniform_int_distribution<unsigned > dist(0, allIndexSwaps.size() - 1);

    for (unsigned i = 0; i < 10; ++i) {
        auto randomIndexSwap = allIndexSwaps[dist(engine)];
        spdlog::debug("Picked permutation: {}",ToString::stdvectorIntToString(randomIndexSwap));
        Eigen::Map<Eigen::VectorXi> temp(randomIndexSwap.data(), randomIndexSwap.size());
        Eigen::PermutationMatrix<Eigen::Dynamic> randomIndexSwapPerm(temp);

        auto Acopy = A;
        Acopy.electrons().permute(randomIndexSwapPerm);

        std::vector<Eigen::VectorXi> permIndices(2, Eigen::VectorXi(B.electrons().numberOfEntities()));
        permIndices[0] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane
        permIndices[1] << 1, 2, 0, 4, 5, 3; // 120° rotation around z

        permIndices[0] = randomIndexSwapPerm * permIndices[0];
        permIndices[1] = randomIndexSwapPerm * permIndices[1];

        MolPermList expectedPerms{
                {{0,1, 2,4,3, 5,7,6}, std::vector<int>(permIndices[0].data(), permIndices[0].data()+permIndices[0].size())},
                {{0,1, 3,4,2, 6,7,5}, std::vector<int>(permIndices[1].data(), permIndices[1].data()+permIndices[1].size())}
        };

        // sort the permutations
        std::sort(std::begin(expectedPerms), std::end(expectedPerms),
                [](const std::pair<std::vector<int>, std::vector<int>>& a,
                        const std::pair<std::vector<int>, std::vector<int>>& b) {
            assert(a.first.size() == b.first.size());
            assert(a.second.size() == b.second.size());

            for (size_t i = 0; i < a.first.size(); ++i)
                if (a.first[i] != b.first[i])
                    return a.first[i] < b.first[i];

            for (size_t i = 0; i < a.second.size(); ++i)
                if (a.second[i] != b.second[i])
                    return a.second[i] < b.second[i];

            return a.first[0] < b.first[0];
        });

        routine(Acopy, B, expectedPerms, distanceTolerance, soapThreshold);
    }
}

TEST_F(AEnvironmentBasedMetricTest, EthaneDoublyIonicAntiMinimal) {
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

    MolPermList expectedPerms{
            {{0,1, 2,4,3, 5,7,6},{0, 2, 1, 3, 5, 4}},// reflection along H2-C0-C1-H5 plane
            {{0,1, 3,4,2, 6,7,5},{1, 2, 0, 4, 5, 3}} // 120° rotation around z
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, EthaneDoublyIonicAnti) {
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

    MolPermList expectedPerms{
            {{0,1, 2,4,3, 5,7,6},{0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16}},// reflection along H2-C0-C1-H5 plane
            {{0,1, 3,4,2, 6,7,5},{0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15}} // 120° rotation around z
    };

    routine(A,B,expectedPerms,distanceTolerance, soapThreshold);
}

TEST_F(AEnvironmentBasedMetricTest, DISABLED_Trans13ButadieneRealMaxima) { // takes too long and fails on GNU/Intel
    auto A = TestMolecules::trans13Butadiene::realA;
    auto B = TestMolecules::trans13Butadiene::realB;

    General::settings.zeta = 3;
    Radial::settings.nmax = 2;
    Angular::settings.lmax = 2;
    Cutoff::settings.radius = 6.0;
    Cutoff::settings.width = 0.0;

    // (permutations are not checked)

    MolPermList expectedPerms{
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}}, // checked via inPsights
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,2,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,1,28,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,28,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,2,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,28,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,2,29}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,2,9}},
    {{3,2,1,0, 9,8,7,6,5,4},{0,18,28,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,2,9}}
    };

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPerms, distanceTolerance, 0.99, true);
}