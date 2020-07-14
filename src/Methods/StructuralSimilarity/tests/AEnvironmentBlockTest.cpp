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


TEST(AEnvironmentBlockTest, FilteringOfCorrectPermutationsAlchemical) {
    auto A = TestMolecules::H4::linear::ionicA;
    auto B = TestMolecules::H4::linear::ionicAreflectedReorderedNumbering; // check with B in chemical mode
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    General::settings.mode = General::Mode::alchemical;
    ParticleKit::create(A);

    std::deque<BestMatch::SOAPSimilarity::GrowingPerm> possiblePermutations = {
            BestMatch::SOAPSimilarity::GrowingPerm({},{{0,1},{2,3}}),
            BestMatch::SOAPSimilarity::GrowingPerm({},{{2,1},{0,3}})
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

    std::deque<BestMatch::SOAPSimilarity::GrowingPerm> possiblePermutations = {
            // second swap is wrong
            BestMatch::SOAPSimilarity::GrowingPerm({},{{0,1},{3,2},{2,3}}),
            // second swap is wrong
            BestMatch::SOAPSimilarity::GrowingPerm({},{{2,1},{3,2},{0,3}})
    };

    auto block = EnvironmentBlock(possiblePermutations, A, B);

    // referenceIndices are ordered
    ASSERT_THAT(block.referenceIndices_, ElementsAre(1,2,3));

    auto filteredPerms = block.filterPermutations(0.1);
    ASSERT_TRUE(filteredPerms.empty());
}
