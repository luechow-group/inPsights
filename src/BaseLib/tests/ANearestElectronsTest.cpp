/* Copyright (C) 2019 Leonard Reuter.
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
#include <TestMolecules.h>
#include <NearestElectrons.h>
#include <Metrics.h>

using namespace testing;

class ANearestElectronsTest : public Test {
public:
    const ElectronsVector &electrons = TestMolecules::BH3::ionic.electrons();
    const AtomsVector &nuclei = TestMolecules::BH3::ionic.atoms();

    void SetUp() override {}
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[0]);
    ASSERT_THAT(indices, ElementsAre(0,4));

    indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[1]);
    ASSERT_THAT(indices,ElementsAre());
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei);
    ASSERT_THAT(indices, ElementsAre(0,4));
};

TEST_F(ANearestElectronsTest, GetElectronsByPositionInverted) {
    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.7);
    auto indices = NearestElectrons::getNearestElectronsIndices(electrons, position, 2);
    ASSERT_THAT(indices, ElementsAre(2,6));
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::function<double(const Eigen::Vector3d &,const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;

    auto indices = NearestElectrons::getNearestElectronsIndices(electrons, nuclei, {nuclei[0].position()}, 2, true,
                                                                20.0, distanceFunction);
    ASSERT_THAT(indices, ElementsAre(6,7));
};

TEST_F(ANearestElectronsTest, InvertedIndices) {
    std::list<long> indices = {0,2,4,6};
    ASSERT_THAT(NearestElectrons::invertedIndices(indices, 7), ElementsAre(1,3,5));
    ASSERT_THAT(NearestElectrons::invertedIndices(indices, 8), ElementsAre(1,3,5,7));

    ASSERT_THAT(NearestElectrons::invertedIndices(std::list<long>{0,1,2}, 3), ElementsAre());

    EXPECT_DEATH(NearestElectrons::invertedIndices(std::list<long>{}, 1),"");
    EXPECT_DEATH(NearestElectrons::invertedIndices(std::list<long>{0}, 0),"");
    EXPECT_DEATH(NearestElectrons::invertedIndices(std::list<long>{0,1,2}, -1),"");
    EXPECT_DEATH(NearestElectrons::invertedIndices(std::list<long>{0,-1,2}, 3),"");
};

TEST_F(ANearestElectronsTest, PickElements) {
    ElectronsVector reference;
    reference.append(electrons[2]); // core electrons are included
    reference.append(electrons[6]);

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.75);
    auto indices = NearestElectrons::getNearestElectronsIndices(electrons, position, 2);

    ASSERT_THAT(indices, ElementsAre(2,6));
    ASSERT_EQ(reference, electrons[indices]);
};

TEST_F(ANearestElectronsTest, GetRelevantIndices) {
    NearestElectrons::settings.distanceMode = "minimum";
    NearestElectrons::settings.maximalCount = 2;
    NearestElectrons::settings.positions = {nuclei[3].position()};
    NearestElectrons::settings.valenceOnly = true;

    auto indices = NearestElectrons::getRelevantIndices(electrons);

    ASSERT_THAT(indices, ElementsAre(3,7));
};
