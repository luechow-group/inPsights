/* Copyright (C) 2020 Michael Heuer.
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
#include <EnvironmentBlock.h>
#include <TestMolecules.h>
#include <ParticleKit.h>

#include <BestMatchSimilarity.h>
#include <SOAPSettings.h>

using namespace testing;
using namespace SOAP;


TEST(AEnvironmentBlockTest, CorrectFilteringAlchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto C = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    // Dependent indices blocks of AC
    // (3 0)
    // (0 1) (2 3)
    // (1 2)

    // pairs of dependent indices
    std::deque<std::pair<Eigen::Index,Eigen::Index>> pairs {
        std::make_pair<Eigen::Index, Eigen::Index>(2,3),
        std::make_pair<Eigen::Index, Eigen::Index>(0,1)};

    EnvironmentBlock block1(pairs, A.electrons(), C.electrons());

    // referenceIndices are ordered
    ASSERT_THAT(block1.referenceIndices_, ElementsAre(1,3));

    // check perm
    auto unfilteredPerms = block1.permutedPermuteeIndicesCollection_;
    ASSERT_THAT(unfilteredPerms[0], ElementsAre(0,2));
    ASSERT_THAT(unfilteredPerms[1], ElementsAre(2,0));

    // filter perm
    auto filteredPerms = block1.filterPermutations(0.1);
    ASSERT_THAT(filteredPerms[0], ElementsAre(0,2));
    ASSERT_THAT(filteredPerms[1], ElementsAre(2,0));
}


TEST(AEnvironmentBlockTest, WrongWithAdditionalPair) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto C = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    // Dependent indices blocks of AC
    // (3 0)
    // (0 1) (2 3)
    // (1 2)

    // pairs of dependent indices
    std::deque<std::pair<Eigen::Index,Eigen::Index>> pairs {
            std::make_pair<Eigen::Index, Eigen::Index>(0,1),
            std::make_pair<Eigen::Index, Eigen::Index>(2,3),
            // wrong index pair:
            std::make_pair<Eigen::Index, Eigen::Index>(1,2),
                    };

    EnvironmentBlock block1(pairs, A.electrons(), C.electrons());

    // referenceIndices are ordered
    ASSERT_THAT(block1.referenceIndices_, ElementsAre(1,2,3));

    // check perm
    auto unfilteredPerms = block1.permutedPermuteeIndicesCollection_;
    ASSERT_EQ(unfilteredPerms.size(), 6);
    ASSERT_THAT(unfilteredPerms[0], ElementsAre(0, 1, 2));
    ASSERT_THAT(unfilteredPerms[1], ElementsAre(0, 2, 1));
    ASSERT_THAT(unfilteredPerms[2], ElementsAre(1, 0, 2));
    ASSERT_THAT(unfilteredPerms[3], ElementsAre(1, 2, 0));
    ASSERT_THAT(unfilteredPerms[4], ElementsAre(2, 0, 1));
    ASSERT_THAT(unfilteredPerms[5], ElementsAre(2, 1, 0));

    auto filteredPerms = block1.filterPermutations(0.1);
    ASSERT_EQ(filteredPerms.size(), 2);
    ASSERT_THAT(filteredPerms[0], ElementsAre(0,1,2)); // TODO understand why they survive
    ASSERT_THAT(filteredPerms[1], ElementsAre(2,1,0));
}

TEST(AEnvironmentBlockTest, UnderstandingEffectOfPermutationOnCovarianceDifference) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto C = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(A);

    // Dependent indices blocks of AC
    // (3 0)
    // (0 1) (2 3)
    // (1 2)

    // In this test, we analyze if (0 1) (2 3) [0] give a different covariance matrix
    // as (2 3) (0 1) [1]
    // or (0 3) (2 1) or (2 1) (0 3)  [2,3]

    auto permutee = A.electrons();
    auto reference = C.electrons();

    std::vector<Eigen::Index> permuteeIndicesOriginal = {0,2};
    std::vector<Eigen::Index> permuteeIndicesOrdered = permuteeIndicesOriginal;
    std::vector<Eigen::Index> permuteeIndicesPermuted = {2,0};

    std::vector<Eigen::Index> referenceIndicesOriginal = {1,3};
    std::vector<Eigen::Index> referenceIndicesOrdered = referenceIndicesOriginal;
    std::vector<Eigen::Index> referenceIndicesPermuted = {3,1};


    spdlog::set_level(spdlog::level::debug);

    double distanceMatrixCovarianceTolerance = 0.01;

    auto orig_orig = DistanceCovariance::conservingQ(
            permuteeIndicesOriginal, permutee,
            referenceIndicesOriginal, reference,
            distanceMatrixCovarianceTolerance);
    spdlog::debug("orig_orig: {}",orig_orig);

    auto perm_orig = DistanceCovariance::conservingQ(
            permuteeIndicesPermuted, permutee,
            referenceIndicesOriginal, reference,
            distanceMatrixCovarianceTolerance);
    spdlog::debug("perm_orig: {}",perm_orig);

    auto orig_perm = DistanceCovariance::conservingQ(
            permuteeIndicesOriginal, permutee,
            referenceIndicesPermuted, reference,
            distanceMatrixCovarianceTolerance);
    spdlog::debug("orig_perm: {}",orig_perm);

    auto perm_perm = DistanceCovariance::conservingQ(
            permuteeIndicesPermuted, permutee,
            referenceIndicesPermuted, reference,
            distanceMatrixCovarianceTolerance);
    spdlog::debug("perm_perm: {}",perm_perm);

}
