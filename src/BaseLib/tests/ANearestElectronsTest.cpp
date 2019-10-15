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
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
    const ElectronsVector &electrons = BH3.electrons();
    const ElectronsVector &electrons2 = TestMolecules::BH3::ionicRotated.electrons();
    ElectronsVector electrons3;
    const AtomsVector &nuclei = BH3.atoms();

    void SetUp() override {
        electrons3 =
                ElectronsVector({
                                        electrons2[0],
                                        electrons2[1],
                                        electrons2[2],
                                        electrons2[7],
                                        electrons2[4],
                                        electrons2[5],
                                        electrons2[6],
                                        electrons2[3]
                                });
    }

    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                               const Eigen::Vector3d &position,
                               const long &count) {
        std::function<double(const Eigen::Vector3d &,
                             const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;
        return NearestElectrons::getNearestElectronsIndices(electrons, nuclei, std::vector<Eigen::Vector3d>({position}),
                                                            count, 100.0,
                                                            distanceFunction, true);
    };
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[0]);
    std::list<long> reference = std::list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[1]);
    reference = std::list<long>({});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei);
    std::list<long> reference({0, 4});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByPosition) {
    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.7);
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electrons, position, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::list<long> indices = getNearestElectronsIndices(electrons, nuclei,
                                                         nuclei[0].position(), 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, PickElements) {
    ElectronsVector reference;
    reference.append(electrons[6]);
    reference.append(electrons[2]);

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.5);

    std::list<long> indices = getNearestElectronsIndices(electrons, nuclei, position, 2);

    ASSERT_EQ(reference, electrons[indices]);
};