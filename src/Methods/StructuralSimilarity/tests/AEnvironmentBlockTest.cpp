// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <EnvironmentBlock.h>
#include <TestMolecules.h>
#include <ParticleKit.h>
#include <EnvironmentBasedMetric.h>
#include <SOAPSettings.h>

using namespace testing;
using namespace SOAP;


TEST(AEnvironmentBlockTest, FilteringOfCorrectPermutationsAlchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    std::deque<Metrics::Similarity::EnvironmentBased::GrowingPerm> possiblePermutations = {
            Metrics::Similarity::EnvironmentBased::GrowingPerm({}, {{0, 1}, {2, 3}}),
            Metrics::Similarity::EnvironmentBased::GrowingPerm({}, {{2, 1}, {0, 3}})
    };

    auto block = EnvironmentBlock(possiblePermutations, A, B);

    // referenceIndices are ordered
    ASSERT_THAT(block.referenceIndices_, ElementsAre(1,3));

    auto filteredPerms = block.filterPermutations(0.1);

    // check perm
    ASSERT_THAT(filteredPerms[0], ElementsAre(0,2));
    ASSERT_THAT(filteredPerms[1], ElementsAre(2,0));
}

TEST(AEnvironmentBlockTest, FilteringOfWrongPermutationsAlchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode

    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    std::deque<Metrics::Similarity::EnvironmentBased::GrowingPerm> possiblePermutations = {
            // second swap is wrong
            Metrics::Similarity::EnvironmentBased::GrowingPerm({}, {{0, 1}, {3, 2}, {2, 3}}),
            // second swap is wrong
            Metrics::Similarity::EnvironmentBased::GrowingPerm({}, {{2, 1}, {3, 2}, {0, 3}})
    };

    auto block = EnvironmentBlock(possiblePermutations, A, B);

    // referenceIndices are ordered
    ASSERT_THAT(block.referenceIndices_, ElementsAre(1,2,3));

    auto filteredPerms = block.filterPermutations(0.1);
    ASSERT_TRUE(filteredPerms.empty());
}
