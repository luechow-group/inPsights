/* Copyright (C) 2019, 2020 Michael Heuer.
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
#include <Motifs.h>
#include <TestMolecules.h>

using namespace testing;

class AMotifsTest : public ::testing::Test {
public:
    Eigen::MatrixXb LiH_filteredSpinCorr, H2_filteredSpinCorr, BH3_filteredSpinCorr;

    Motifs LiH, H2covalent, H2ionicRight, BH3covalent;

    void SetUp() override {

        MolecularGeometry molLiH = {
                AtomsVector({
                    {Element::Li, {0, 0, 1}},
                    {Element::H,  {0, 0, -1}}}),
                ElectronsVector({
                    {Spin::alpha, {0, 0, 1}},
                    {Spin::beta,  {0, 0, 1}},
                    {Spin::alpha, {0, 0, -1}},
                    {Spin::beta,  {0, 0, -0.5}},
                    })};

        LiH_filteredSpinCorr = Eigen::MatrixXb(4, 4);
        LiH_filteredSpinCorr <<
        0, 1, 0, 0, \
        1, 0, 0, 0, \
        0, 0, 0, 1, \
        0, 0, 1, 0;
        LiH = Motifs(LiH_filteredSpinCorr, molLiH);

        H2_filteredSpinCorr = Eigen::MatrixXb(2, 2);
        H2_filteredSpinCorr <<
        0, 1, \
        1, 0;
        H2covalent = Motifs(H2_filteredSpinCorr, TestMolecules::H2::ElectronsInCores::normal);
        H2ionicRight = Motifs(H2_filteredSpinCorr, TestMolecules::H2::ElectronsInCores::ionicRight);


        BH3_filteredSpinCorr = Eigen::MatrixXb(8, 8);
        BH3_filteredSpinCorr <<
        0, 0, 0, 0, 1, 0, 0, 0, \
        0, 0, 0, 0, 0, 1, 0, 0, \
        0, 0, 0, 0, 0, 0, 1, 0, \
        0, 0, 0, 0, 0, 0, 0, 1, \
        1, 0, 0, 0, 0, 0, 0, 0, \
        0, 1, 0, 0, 0, 0, 0, 0, \
        0, 0, 1, 0, 0, 0, 0, 0, \
        0, 0, 0, 1, 0, 0, 0, 0;
        BH3covalent = Motifs(BH3_filteredSpinCorr, TestMolecules::BH3::normal);
    };
};

TEST_F(AMotifsTest, LiH) {

    ASSERT_THAT(LiH.motifVector_[0].electronIndices(), ElementsAre(0, 1));
    ASSERT_THAT(LiH.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(LiH.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(LiH.motifVector_[1].electronIndices(), ElementsAre());
    ASSERT_THAT(LiH.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(LiH.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(LiH.motifVector_[2].electronIndices(), ElementsAre(2, 3));
    ASSERT_THAT(LiH.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(LiH.motifVector_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, H2covalent) {

    ASSERT_THAT(H2covalent.motifVector_[0].electronIndices(), ElementsAre());
    ASSERT_THAT(H2covalent.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(H2covalent.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(H2covalent.motifVector_[1].electronIndices(), ElementsAre());
    ASSERT_THAT(H2covalent.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(H2covalent.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(H2covalent.motifVector_[2].electronIndices(), ElementsAre(0, 1));
    ASSERT_THAT(H2covalent.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(H2covalent.motifVector_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, H2ionicRight) {

    ASSERT_THAT(H2ionicRight.motifVector_[0].electronIndices(), ElementsAre());
    ASSERT_THAT(H2ionicRight.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(H2ionicRight.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(H2ionicRight.motifVector_[1].electronIndices(), ElementsAre());
    ASSERT_THAT(H2ionicRight.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(H2ionicRight.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(H2ionicRight.motifVector_[2].electronIndices(), ElementsAre(0, 1));
    ASSERT_THAT(H2ionicRight.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(H2ionicRight.motifVector_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, BH3covalent) {
    ASSERT_THAT(BH3covalent.motifVector_[0].electronIndices(), ElementsAre(0, 4));
    ASSERT_THAT(BH3covalent.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(BH3covalent.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifVector_[1].electronIndices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(BH3covalent.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifVector_[2].electronIndices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifVector_[2].atomIndices(), ElementsAre(2));
    ASSERT_EQ(BH3covalent.motifVector_[2].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifVector_[3].electronIndices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifVector_[3].atomIndices(), ElementsAre(3));
    ASSERT_EQ(BH3covalent.motifVector_[3].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifVector_[4].electronIndices(), ElementsAre(1, 5));
    ASSERT_THAT(BH3covalent.motifVector_[4].atomIndices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifVector_[4].type(), MotifType::Valence);

    ASSERT_THAT(BH3covalent.motifVector_[5].electronIndices(), ElementsAre(2, 6));
    ASSERT_THAT(BH3covalent.motifVector_[5].atomIndices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifVector_[5].type(), MotifType::Valence);

    ASSERT_THAT(BH3covalent.motifVector_[6].electronIndices(), ElementsAre(3, 7));
    ASSERT_THAT(BH3covalent.motifVector_[6].atomIndices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifVector_[6].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, YAML) {
    auto node = YAML::convert<Motif>::encode(LiH.motifVector_[0]);

    Motif decodedMotif;
    YAML::convert<Motif>::decode(node, decodedMotif);

    ASSERT_THAT(decodedMotif.electronIndices(), ElementsAre(0, 1));
    ASSERT_THAT(decodedMotif.atomIndices(), ElementsAre(0));
    ASSERT_EQ(decodedMotif.type(), MotifType::Core);
}
