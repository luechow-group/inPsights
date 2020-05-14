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
#include <limits>
#include <random>
#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <Metrics.h>
#include <Hungarian.h>
#include <Combinatorics.h>
#include <Enumerate.h>

using namespace testing;
using namespace SOAP;

class ABestMatchSimilarityTest : public ::testing::Test {
public:
    double distanceTolerance, soapThreshold, shakeSoapThreshold, numericalPrecisionEpsilon;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        //spdlog::set_level(spdlog::level::debug);

        distanceTolerance = 0.1;
        soapThreshold = 1.0;
        shakeSoapThreshold = 0.90;
        numericalPrecisionEpsilon = SOAP::General::settings.comparisonEpsilon.get();

        Radial::settings.nmax = 2;
        Radial::settings.sigmaAtom = 2.0;
        Angular::settings.lmax = 2;
        Cutoff::settings.radius = 6.0;
        Cutoff::settings.width = 4.0;
        General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    };

    void routine(const MolecularGeometry &A, const MolecularGeometry &B,
                 const std::vector<Eigen::VectorXi> &expectedPermutationIndices,
                 double distTolerance, double soapThresh, bool greaterThan = false) {

        ParticleKit::create(A);
        ASSERT_TRUE(ParticleKit::isSubsetQ(A));
        ASSERT_TRUE(ParticleKit::isSubsetQ(B));

        auto specA = MolecularSpectrum(A);
        auto specB = MolecularSpectrum(B);

        auto results = BestMatch::SOAPSimilarity::getBestMatchResults(
                specA, specB,
                distTolerance, soapThresh, numericalPrecisionEpsilon);

        std::sort(results.begin(), results.end());

        spdlog::debug("Expected:");
        for (auto[i, expected] : enumerate(expectedPermutationIndices))
            spdlog::debug("{}: metric {} {}, permutation = {} (lab system)", i, greaterThan ? ">" : "=", soapThresh,
                          ToString::vectorXiToString(expected));


        unsigned removedCounter = 0;
        for(auto it = results.begin(); it != results.end(); it++) {
            if(it->metric < (soapThresh-numericalPrecisionEpsilon)) {
                results.erase(it--);
                removedCounter++;
            }
        }
        spdlog::debug("Removed {} result(s) below threshold.", removedCounter);

        spdlog::debug("Got:");
        for (auto[i, result] : enumerate(results)) {
            spdlog::debug("{}: metric = {:01.16f}, permutation = {} (lab system)", i, result.metric,
                          ToString::vectorXiToString(result.permutation.indices()));
        }

        size_t i = 0;
        while (i < results.size() || i < expectedPermutationIndices.size()) {

            if (i >= results.size()) {
                spdlog::debug(" {}: missing", i);
                i++;
                continue;
            }

            if (i < expectedPermutationIndices.size()) {
                if (!greaterThan) ASSERT_NEAR(results[i].metric, soapThresh, numericalPrecisionEpsilon);
                else ASSERT_GT(results[i].metric, soapThresh);

                ASSERT_EQ(results[i].permutation.indices(), expectedPermutationIndices[i]);
            } else {
                spdlog::debug("Found unexpected result {}", i);
            };
            i++;

        }
        ASSERT_EQ(results.size(), expectedPermutationIndices.size());
    }
};

TEST_F(ABestMatchSimilarityTest, PrermuteEnvironmentsToLabSystem) {
    General::settings.mode = General::Mode::chemical;

    auto B = TestMolecules::BH3::ionicMinimal;
    auto A = TestMolecules::BH3::ionicMinimalRotatedPermuted;
    ParticleKit::create(A);

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto environmentalSimilarities = BestMatch::SOAPSimilarity::calculateEnvironmentSimilarityMatrix(specA, specB);

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

TEST_F(ABestMatchSimilarityTest, GrowingPerm) {

    auto growingPerm = BestMatch::SOAPSimilarity::GrowingPerm({0,1,2},{});

    ASSERT_TRUE(growingPerm.add({0,1}));
    ASSERT_THAT(growingPerm.remainingPermuteeIndices_,ElementsAre(1,2));

    ASSERT_EQ(growingPerm.chainOfSwaps_[0].first, 0);
    ASSERT_EQ(growingPerm.chainOfSwaps_[0].second, 1);

    ASSERT_FALSE(growingPerm.add({0,2}));
}

TEST_F(ABestMatchSimilarityTest, GroupDependentMatches1) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,0,1,0,\
    0,1,0,1,\
    1,0,1,0,\
    0,1,0,1;

    auto matches = BestMatch::SOAPSimilarity::findEnvironmentMatches(mat, 1.0, 0.001);
    ASSERT_THAT(matches[0].permuteeEnvsIndices, ElementsAre(0,2));
    ASSERT_EQ(matches[0].referenceEnvIndex, 0);
    ASSERT_THAT(matches[1].permuteeEnvsIndices, ElementsAre(1,3));
    ASSERT_EQ(matches[1].referenceEnvIndex, 1);
    ASSERT_THAT(matches[2].permuteeEnvsIndices, ElementsAre(0,2));
    ASSERT_EQ(matches[2].referenceEnvIndex, 2);
    ASSERT_THAT(matches[3].permuteeEnvsIndices, ElementsAre(1,3));
    ASSERT_EQ(matches[3].referenceEnvIndex, 3);

    auto dependentMatches = BestMatch::SOAPSimilarity::groupDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeEnvsIndices, ElementsAre(0,2));
    ASSERT_EQ(dependentMatches[0][0].referenceEnvIndex, 0);
    ASSERT_THAT(dependentMatches[0][1].permuteeEnvsIndices, ElementsAre(0,2));
    ASSERT_EQ(dependentMatches[0][1].referenceEnvIndex, 2);
    ASSERT_THAT(dependentMatches[1][0].permuteeEnvsIndices, ElementsAre(1,3));
    ASSERT_EQ(dependentMatches[1][0].referenceEnvIndex, 1);
    ASSERT_THAT(dependentMatches[1][1].permuteeEnvsIndices, ElementsAre(1,3));
    ASSERT_EQ(dependentMatches[1][1].referenceEnvIndex, 3);
}

TEST_F(ABestMatchSimilarityTest, GroupDependentMatches2) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,1,0,1,\
    1,1,0,0,\
    0,0,1,0,\
    1,0,0,1;

    auto matches = BestMatch::SOAPSimilarity::findEnvironmentMatches(mat, 1.0, 1e-8);
    ASSERT_THAT(matches[0].permuteeEnvsIndices, ElementsAre(0,1,3));
    ASSERT_EQ(matches[0].referenceEnvIndex, 0);
    ASSERT_THAT(matches[1].permuteeEnvsIndices, ElementsAre(0,1));
    ASSERT_EQ(matches[1].referenceEnvIndex, 1);
    ASSERT_THAT(matches[2].permuteeEnvsIndices, ElementsAre(2));
    ASSERT_EQ(matches[2].referenceEnvIndex, 2);
    ASSERT_THAT(matches[3].permuteeEnvsIndices, ElementsAre(0,3));
    ASSERT_EQ(matches[3].referenceEnvIndex, 3);

    auto dependentMatches = BestMatch::SOAPSimilarity::groupDependentMatches(matches);
    ASSERT_THAT(dependentMatches[0][0].permuteeEnvsIndices, ElementsAre(0,1,3));
    ASSERT_EQ(dependentMatches[0][0].referenceEnvIndex, 0);
    ASSERT_THAT(dependentMatches[0][1].permuteeEnvsIndices, ElementsAre(0,1));
    ASSERT_EQ(dependentMatches[0][1].referenceEnvIndex, 1);
    ASSERT_THAT(dependentMatches[0][2].permuteeEnvsIndices, ElementsAre(0,3));
    ASSERT_EQ(dependentMatches[0][2].referenceEnvIndex, 3);

    ASSERT_THAT(dependentMatches[1][0].permuteeEnvsIndices, ElementsAre(2));
    ASSERT_EQ(dependentMatches[1][0].referenceEnvIndex, 2);

    auto possiblePermutations = BestMatch::SOAPSimilarity::findPossiblePermutations(dependentMatches[0]);

    using Pair = std::pair<Eigen::Index,Eigen::Index>;
    ASSERT_EQ(possiblePermutations.size(), 3);
    ASSERT_TRUE(possiblePermutations[0].remainingPermuteeIndices_.empty());
    ASSERT_THAT(possiblePermutations[0].chainOfSwaps_, ElementsAre(Pair(3,0),Pair(1,1),Pair(0,3)));

    ASSERT_TRUE(possiblePermutations[1].remainingPermuteeIndices_.empty());
    ASSERT_THAT(possiblePermutations[1].chainOfSwaps_, ElementsAre(Pair(1,0),Pair(0,1),Pair(3,3)));

    ASSERT_TRUE(possiblePermutations[2].remainingPermuteeIndices_.empty());
    ASSERT_THAT(possiblePermutations[2].chainOfSwaps_, ElementsAre(Pair(0,0),Pair(1,1),Pair(3,3)));

}



TEST_F(ABestMatchSimilarityTest, ListOfDependentIndicesCase1) {
    Eigen::MatrixXd mat(4,4);
    mat <<
    1,0,1,0,\
    0,1,0,1,\
    1,0,1,0,\
    0,1,0,1;

    auto bestMatch = Hungarian<double>::findMatching(mat, Matchtype::MAX);

    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(mat.rows());
    expectedPerm.setIdentity();  // rows are permuted
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto bestMatchPermutedEnvironments = bestMatch*mat;
    auto result = BestMatch::SOAPSimilarity::findEquivalentEnvironments(bestMatchPermutedEnvironments, bestMatch, 1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_EQ(result[0][0].first, 0);
    ASSERT_EQ(result[0][0].second, 0);
    ASSERT_EQ(result[0][1].first, 2);
    ASSERT_EQ(result[0][1].second, 2);

    ASSERT_EQ(result[1][0].first, 1);
    ASSERT_EQ(result[1][0].second, 1);
    ASSERT_EQ(result[1][1].first, 3);
    ASSERT_EQ(result[1][1].second, 3);
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndicesCase2) {
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

    auto bestMatchPermutedEnvironments = bestMatch*mat;
    auto result = BestMatch::SOAPSimilarity::findEquivalentEnvironments(bestMatchPermutedEnvironments, bestMatch, 1.0);

    ASSERT_EQ(result.size(),2);
    ASSERT_EQ(result[0][0].first, 0);
    ASSERT_EQ(result[0][0].second, 0);
    ASSERT_EQ(result[0][1].first, 1);
    ASSERT_EQ(result[0][1].second, 1);

    ASSERT_EQ(result[1][0].first, 2);
    ASSERT_EQ(result[1][0].second, 2);
    ASSERT_EQ(result[1][1].first, 3);
    ASSERT_EQ(result[1][1].second, 3);
}

TEST_F(ABestMatchSimilarityTest, ListOfDependentIndicesCase3) {
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
    expectedPermIndices << 0,2,1,3,4,5; // rows are permuted
    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(expectedPermIndices);
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto bestMatchPermutedEnvironments = bestMatch*mat;
    auto result = BestMatch::SOAPSimilarity::findEquivalentEnvironments(bestMatchPermutedEnvironments, bestMatch, 1.0);

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

TEST_F(ABestMatchSimilarityTest, FindEquivalentEnvironmentsCase4) {
    Eigen::MatrixXd mat(6,6);
    mat <<
    0,1,0,0,0,0,\
    0,0,1,0,0,0,\
    1,0,0,0,0,0,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1,\
    0,0,0,1,1,1;

    auto bestMatch = Hungarian<double>::findMatching(mat, Matchtype::MAX);

    Eigen::VectorXi expectedPermIndices(mat.rows());
    expectedPermIndices << 1,2,0,3,4,5; // rows are permuted
    Eigen::PermutationMatrix<Eigen::Dynamic> expectedPerm(expectedPermIndices);
    ASSERT_TRUE(bestMatch.indices().isApprox(expectedPerm.indices()));

    auto bestMatchPermutedEnvironments = bestMatch*mat;
    auto result = BestMatch::SOAPSimilarity::findEquivalentEnvironments(bestMatchPermutedEnvironments, bestMatch, 1.0);

    ASSERT_EQ(result.size(),4);

    ASSERT_EQ(result[0].size(),1);
    ASSERT_EQ(result[0][0].first, 2);
    ASSERT_EQ(result[0][0].second,0);

    ASSERT_EQ(result[1].size(),1);
    ASSERT_EQ(result[1][0].first, 0);
    ASSERT_EQ(result[1][0].second,1);

    ASSERT_EQ(result[2].size(),1);
    ASSERT_EQ(result[2][0].first, 1);
    ASSERT_EQ(result[2][0].second,2);

    ASSERT_EQ(result[3].size(),3);
    ASSERT_EQ(result[3][0].first, 3);
    ASSERT_EQ(result[3][0].second,3);
    ASSERT_EQ(result[3][1].first, 4);
    ASSERT_EQ(result[3][1].second,4);
    ASSERT_EQ(result[3][2].first, 5);
    ASSERT_EQ(result[3][2].second,5);
}


TEST_F(ABestMatchSimilarityTest, H4linear_alchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto C = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 1,2,3,0;
    expectedPermIndices[1] << 3,2,1,0;
    routine(A, C, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_ionic_reflected) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflected;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,2,3;
    expectedPermIndices[1] << 2,1,0,3;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPermIndices = {expectedPermIndices[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_ionic_reflected_alpha_permuted) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedAlphaPermuted;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 1,0,2,3;
    expectedPermIndices[1] << 2,0,1,3;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPermIndices = {expectedPermIndices[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_ionic_reflected_beta_permuted) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedBetaPermuted;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,3,2;
    expectedPermIndices[1] << 3,1,0,2;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPermIndices = {expectedPermIndices[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_ionic_reflected_reordered_numbering) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 1,2,3,0;
    expectedPermIndices[1] << 3,2,1,0;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    // both perms are not conserving in chemical mode
    expectedPermIndices = {};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_not_in_particle_kit) {
    auto B = TestMolecules::H4::linear::ionicA;
    auto A = TestMolecules::H4::linear::ionicNotinParticleKitSystem;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,2,3;
    expectedPermIndices[1] << 2,1,0,3;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    // second perm is not conserving in chemical mode
    expectedPermIndices = {expectedPermIndices[1]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_real_maxima1) {
    auto A = TestMolecules::H4::linear::ionicRealMax1Var1;
    auto B = TestMolecules::H4::linear::ionicRealMax1Var2;

    std::vector<Eigen::VectorXi> expectedPermIndices(2, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,2,1,3;
    expectedPermIndices[1] << 3,2,1,0;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    General::settings.mode = General::Mode::chemical;
    auto envSimMat = BestMatch::SOAPSimilarity::calculateEnvironmentSimilarityMatrix(MolecularSpectrum(A),
                                                                                     MolecularSpectrum(B));
    ASSERT_GT(envSimMat.minCoeff(), 0.0);
    ASSERT_LT(envSimMat.maxCoeff(), 1.0); // no equivalent environments exist in chemical mode.
    expectedPermIndices = {};
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4linear_real_maxima2) {
    auto A = TestMolecules::H4::linear::ionicRealMax2Var1;
    auto B = TestMolecules::H4::linear::ionicRealMax2Var2;

    std::vector<Eigen::VectorXi> expectedPermIndices(4, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 1,0,3,2;
    expectedPermIndices[1] << 1,3,0,2;
    expectedPermIndices[2] << 2,0,3,1;
    expectedPermIndices[3] << 2,3,0,1;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);

    expectedPermIndices = {expectedPermIndices[0]};
    General::settings.mode = General::Mode::chemical;
    routine(A, B, expectedPermIndices, distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4ring) {
    auto B = TestMolecules::H4::ring::fourAlpha;
    auto A = B;

    std::vector<Eigen::VectorXi> permsIndices(8,Eigen::VectorXi(B.electrons().numberOfEntities()));
    permsIndices[0] << 0,1,2,3;
    permsIndices[1] << 0,1,3,2;
    permsIndices[2] << 1,0,2,3;
    permsIndices[3] << 1,0,3,2;
    permsIndices[4] << 2,3,0,1;
    permsIndices[5] << 2,3,1,0;
    permsIndices[6] << 3,2,0,1;
    permsIndices[7] << 3,2,1,0;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A,B,permsIndices,distanceTolerance, soapThreshold);

    General::settings.mode = General::Mode::chemical;
    routine(A,B,permsIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, H4ring_Shaked) {
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

        std::vector<Eigen::VectorXi> permsIndices(8, Eigen::VectorXi(B.electrons().numberOfEntities()));
        permsIndices[0] << 0,1,2,3;
        permsIndices[1] << 0,1,3,2;
        permsIndices[2] << 1,0,2,3;
        permsIndices[3] << 1,0,3,2;
        permsIndices[4] << 2,3,0,1;
        permsIndices[5] << 2,3,1,0;
        permsIndices[6] << 3,2,0,1;
        permsIndices[7] << 3,2,1,0;

        General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
        General::settings.mode = General::Mode::alchemical;
        routine(A, B, permsIndices, distanceTolerance, shakeSoapThreshold, true);

        General::settings.mode = General::Mode::chemical;
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
    permIndices[0] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane
    permIndices[1] << 1, 2, 0, 4, 5, 3; // 120° rotation around z

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}


TEST_F(ABestMatchSimilarityTest, EthaneSinglyIonicMinimal_shaked_alchemical) {
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
    permIndices[0] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane
    permIndices[1] << 1, 2, 0, 4, 5, 3; // 120° rotation around z


    //auto randomSeed = static_cast<unsigned long>(123);
    //std::cout << "random seed: " << randomSeed << std::endl;

    auto rng = std::default_random_engine(123);
    std::cout<< std::endl << A.electrons() << std::endl;
    A.electrons().positionsVector().shake(distanceTolerance / 10.0, rng);
    std::cout << A.electrons() << std::endl;

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    routine(A, B, permIndices, distanceTolerance, shakeSoapThreshold, true);
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
    permIndices[0] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane
    permIndices[1] << 0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15; // 120° rotation around z

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, EthaneSinglyIonic10RandomIndexSwapPermutations) {

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

        // sort the permutations
        std::sort(std::begin(permIndices), std::end(permIndices),
                [](const Eigen::VectorXi& a,const Eigen::VectorXi& b) {

            assert(a.size() == b.size());
            for (Eigen::Index i = 0; i < a.size(); ++i) {
                if (a[i] != b[i])
                    return a[i] < b[i];
            }
            return a[0] < b[0];
        });

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
    permIndices[0] << 0, 2, 1, 3, 5, 4; // reflection along H2-C0-C1-H5 plane
    permIndices[1] << 1, 2, 0, 4, 5, 3; // 120° rotation around z

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
    permIndices[0] << 0,1,2,3,  4,6,5, 7,9,8,  10,11,  12,14,13, 15,17,16; // reflection along H2-C0-C1-H5 plane
    permIndices[1] << 0,1,2,3,  5,6,4, 8,9,7,  10,11,  13,14,12, 16,17,15; // 120° rotation around z

    routine(A,B,permIndices,distanceTolerance, soapThreshold);
}

TEST_F(ABestMatchSimilarityTest, Trans13ButadieneRealMaxima) {
    spdlog::set_level(spdlog::level::debug);
    auto A = TestMolecules::trans13Butadiene::realA;
    auto B = TestMolecules::trans13Butadiene::realB;

    General::settings.zeta = 3;
    Radial::settings.nmax = 3;
    Angular::settings.lmax = 3;
    Cutoff::settings.radius = 8.0;
    Cutoff::settings.width = 8.0;
    Radial::settings.sigmaAtom = 1.0;

    // (permutations are not checked)
    std::vector<Eigen::VectorXi> expectedPermIndices(64, Eigen::VectorXi(A.electrons().numberOfEntities()));
    expectedPermIndices[0] << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29; // checked via inPsights
    expectedPermIndices[1] << 0,1,2,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,28,29;
    expectedPermIndices[2] << 0,1,2,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,9;
    expectedPermIndices[3] << 0,1,2,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,28,9;
    expectedPermIndices[4] << 0,1,2,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,28,29;
    expectedPermIndices[5] << 0,1,2,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,28,29;
    expectedPermIndices[6] << 0,1,2,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,28,9;
    expectedPermIndices[7] << 0,1,2,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,28,9;
    expectedPermIndices[8] << 0,1,28,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,2,29;
    expectedPermIndices[9] << 0,1,28,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,2,29;
    expectedPermIndices[10] << 0,1,28,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,2,9;
    expectedPermIndices[11] << 0,1,28,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,24,25,26,27,2,9;
    expectedPermIndices[12] << 0,1,28,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,2,29;
    expectedPermIndices[13] << 0,1,28,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,2,29;
    expectedPermIndices[14] << 0,1,28,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,18,19,20,21,22,23,6,25,26,27,2,9;
    expectedPermIndices[15] << 0,1,28,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,18,19,20,21,22,23,6,25,26,27,2,9;
    expectedPermIndices[16] << 0,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,28,29;
    expectedPermIndices[17] << 0,18,2,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,28,29;
    expectedPermIndices[18] << 0,18,2,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,28,9;
    expectedPermIndices[19] << 0,18,2,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,28,9;
    expectedPermIndices[20] << 0,18,2,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,28,29;
    expectedPermIndices[21] << 0,18,2,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,28,29;
    expectedPermIndices[22] << 0,18,2,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,28,9;
    expectedPermIndices[23] << 0,18,2,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,28,9;
    expectedPermIndices[24] << 0,18,28,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,2,29;
    expectedPermIndices[25] << 0,18,28,3,4,5,6,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,2,29;
    expectedPermIndices[26] << 0,18,28,3,4,5,6,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,24,25,26,27,2,9;
    expectedPermIndices[27] << 0,18,28,3,4,5,6,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,24,25,26,27,2,9;
    expectedPermIndices[28] << 0,18,28,3,4,5,24,7,8,9,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,2,29;
    expectedPermIndices[29] << 0,18,28,3,4,5,24,7,8,9,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,2,29;
    expectedPermIndices[30] << 0,18,28,3,4,5,24,7,8,29,10,11,12,13,14,15,16,17,1,19,20,21,22,23,6,25,26,27,2,9;
    expectedPermIndices[31] << 0,18,28,3,4,5,24,7,8,29,10,11,15,13,14,12,16,17,1,19,20,21,22,23,6,25,26,27,2,9;
    expectedPermIndices[32] << 0,1,2,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,28,29;
    expectedPermIndices[33] << 0,1,2,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,28,29;
    expectedPermIndices[34] << 0,1,2,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,28,9;
    expectedPermIndices[35] << 0,1,2,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,28,9;
    expectedPermIndices[36] << 0,1,2,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,28,29;
    expectedPermIndices[37] << 0,1,2,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,28,29;
    expectedPermIndices[38] << 0,1,2,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,28,9;
    expectedPermIndices[39] << 0,1,2,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,28,9;
    expectedPermIndices[40] << 0,1,28,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,2,29;
    expectedPermIndices[41] << 0,1,28,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,2,29;
    expectedPermIndices[42] << 0,1,28,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,24,17,19,27,2,9;
    expectedPermIndices[43] << 0,1,28,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,24,17,19,27,2,9;
    expectedPermIndices[44] << 0,1,28,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,2,29;
    expectedPermIndices[45] << 0,1,28,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,2,29;
    expectedPermIndices[46] << 0,1,28,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,18,26,20,21,22,23,6,17,19,27,2,9;
    expectedPermIndices[47] << 0,1,28,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,18,26,20,21,22,23,6,17,19,27,2,9;
    expectedPermIndices[48] << 0,18,2,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,28,29;
    expectedPermIndices[49] << 0,18,2,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,28,29;
    expectedPermIndices[50] << 0,18,2,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,28,9;
    expectedPermIndices[51] << 0,18,2,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,28,9;
    expectedPermIndices[52] << 0,18,2,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,28,29;
    expectedPermIndices[53] << 0,18,2,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,28,29;
    expectedPermIndices[54] << 0,18,2,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,28,9;
    expectedPermIndices[55] << 0,18,2,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,28,9;
    expectedPermIndices[56] << 0,18,28,11,4,5,6,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,2,29;
    expectedPermIndices[57] << 0,18,28,11,4,5,6,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,2,29;
    expectedPermIndices[58] << 0,18,28,11,4,5,6,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,24,17,19,27,2,9;
    expectedPermIndices[59] << 0,18,28,11,4,5,6,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,24,17,19,27,2,9;
    expectedPermIndices[60] << 0,18,28,11,4,5,24,7,8,9,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,2,29;
    expectedPermIndices[61] << 0,18,28,11,4,5,24,7,8,9,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,2,29;
    expectedPermIndices[62] << 0,18,28,11,4,5,24,7,8,29,14,3,12,13,10,15,16,25,1,26,20,21,22,23,6,17,19,27,2,9;
    expectedPermIndices[63] << 0,18,28,11,4,5,24,7,8,29,14,3,15,13,10,12,16,25,1,26,20,21,22,23,6,17,19,27,2,9;

    General::settings.mode = General::Mode::typeAgnostic;
    routine(A, B, expectedPermIndices, distanceTolerance, 0.99, true);
}