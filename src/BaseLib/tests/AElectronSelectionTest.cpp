// Copyright (C) 2019 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <TestMolecules.h>
#include <ElectronSelection.h>
#include <Metrics.h>

using namespace testing;

class AElectronSelectionTest : public Test {
public:
    const ElectronsVector &electrons = TestMolecules::BH3::ionic.electrons();
    const AtomsVector &nuclei = TestMolecules::BH3::ionic.atoms();

    void SetUp() override {}
};

TEST_F(AElectronSelectionTest, GetNonValenceIndices) {
    std::list<long> indices = ElectronSelection::getNonValenceIndices(electrons, nuclei[0]);
    ASSERT_THAT(indices, ElementsAre(0,4));

    indices = ElectronSelection::getNonValenceIndices(electrons, nuclei[1]);
    ASSERT_THAT(indices,ElementsAre());
};

TEST_F(AElectronSelectionTest, GetNonValenceIndicesAll) {
    std::list<long> indices = ElectronSelection::getNonValenceIndices(electrons, nuclei);
    ASSERT_THAT(indices, ElementsAre(0,4));
};

TEST_F(AElectronSelectionTest, SetTest) {
    std::set<long> indices = {0,4};
    ASSERT_THAT(indices, ElementsAre(0,4));
};

TEST_F(AElectronSelectionTest, GetElectronsByPositionInverted) {
    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.7);
    auto indices = ElectronSelection::getNearestPositionIndices(electrons.positionsVector(), position, 2);
    ASSERT_THAT(indices, ElementsAre(2,6));
};

TEST_F(AElectronSelectionTest, GetValenceByPosition) {
    std::function<double(const Eigen::Vector3d &,const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;

    auto indices = ElectronSelection::getNearestElectronsIndices(electrons, nuclei, {nuclei[0].position()}, 2, true,
                                                                 20.0, distanceFunction);
    ASSERT_THAT(indices, ElementsAre(6,7));
};

TEST_F(AElectronSelectionTest, InvertedIndices) {
    std::list<long> indices = {0,2,4,6};
    ASSERT_THAT(ElectronSelection::invertedIndices(indices, 7), ElementsAre(1, 3, 5));
    ASSERT_THAT(ElectronSelection::invertedIndices(indices, 8), ElementsAre(1, 3, 5, 7));

    ASSERT_THAT(ElectronSelection::invertedIndices(std::list<long>{0, 1, 2}, 3), ElementsAre());

    EXPECT_DEATH(ElectronSelection::invertedIndices(std::list<long>{}, 1), "");
    EXPECT_DEATH(ElectronSelection::invertedIndices(std::list<long>{0}, 0), "");
    EXPECT_DEATH(ElectronSelection::invertedIndices(std::list<long>{0, 1, 2}, -1), "");
    EXPECT_DEATH(ElectronSelection::invertedIndices(std::list<long>{0, -1, 2}, 3), "");
};

TEST_F(AElectronSelectionTest, InvertedIndicesUnordered) {
    std::list<long> indices = {0,4,2,6};
    ASSERT_THAT(ElectronSelection::invertedIndices(indices, 7),
                ElementsAre(1, 3, 5));
    ASSERT_THAT(ElectronSelection::invertedIndices(indices, 8), ElementsAre(1, 3, 5, 7));
};

TEST_F(AElectronSelectionTest, PickElements) {
    ElectronsVector reference;
    reference.append(electrons[2]); // core electrons are included
    reference.append(electrons[6]);

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.75);
    auto indices = ElectronSelection::getNearestPositionIndices(electrons.positionsVector(), position, 2);

    ASSERT_THAT(indices, ElementsAre(2,6));
    ASSERT_EQ(reference, electrons[indices]);
};

TEST_F(AElectronSelectionTest, GetRelevantIndices) {
    ElectronSelection::settings.distanceMode = "minimum";
    ElectronSelection::settings.maximalCount = 2;
    ElectronSelection::settings.positions = {nuclei[3].position()};
    ElectronSelection::settings.valenceOnly = true;

    auto indices = ElectronSelection::getRelevantIndices(electrons, nuclei);

    ASSERT_THAT(indices, ElementsAre(3,7));
};

TEST_F(AElectronSelectionTest, BH3_InvertSelectionValenceOnly) {
    // Expectation:
    // All electrons except the two at the ionic position at nucleus 1

    ElectronSelection::settings.distanceMode = "minimum";
    ElectronSelection::settings.maximalCount = 2;
    ElectronSelection::settings.positions = {nuclei[1].position()};
    ElectronSelection::settings.valenceOnly = false;
    ElectronSelection::settings.invertSelection = true;

    auto indices = ElectronSelection::getRelevantIndices(electrons, nuclei);

    ASSERT_THAT(indices, ElementsAre(0,2,3,4,6,7));
};

TEST_F(AElectronSelectionTest, BH3_InvertSelection) {

    // Expectation:
    // All electrons except the two core electrons in boron and the two at the ionic position at nucleus 1

    ElectronSelection::settings.distanceMode = "minimum";
    ElectronSelection::settings.maximalCount = 2;
    ElectronSelection::settings.positions = {nuclei[1].position()};
    ElectronSelection::settings.valenceOnly = true;
    ElectronSelection::settings.invertSelection = true;

    auto indices = ElectronSelection::getRelevantIndices(electrons, nuclei);

    ASSERT_THAT(indices, ElementsAre(2,3,6,7));
};