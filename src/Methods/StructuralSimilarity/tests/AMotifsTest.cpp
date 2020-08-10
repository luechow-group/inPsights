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
    Eigen::MatrixXb LiH_filteredSpinCorr, H2_filteredSpinCorr, BH3_filteredSpinCorr, Li2_filteredSpinCorr;

    Motifs LiH, H2covalent, H2ionicRight, BH3covalent, Li2;

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

        Li2_filteredSpinCorr = Eigen::MatrixXb(6, 6);
        Li2_filteredSpinCorr <<
        0, 1, 1, 1, 0, 0, \
        1, 0, 1, 1, 0, 0, \
        1, 1, 0, 1, 0, 0, \
        1, 1, 1, 0, 0, 0, \
        0, 0, 0, 0, 0, 1, \
        0, 0, 0, 0, 1, 0;
        Li2 = Motifs(Li2_filteredSpinCorr, TestMolecules::Li2::normal);
    };
};

TEST_F(AMotifsTest, Li2) {
    // The connected first four electrons are split into two core motifs
    ASSERT_THAT(Li2.motifs_[0].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(Li2.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(Li2.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(Li2.motifs_[1].electrons_.indices(), ElementsAre(2, 3));
    ASSERT_THAT(Li2.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(Li2.motifs_[1].type(), MotifType::Core);

    // Valence elctrons are not influenced by this
    ASSERT_THAT(Li2.motifs_[2].electrons_.indices(), ElementsAre(4, 5));
    ASSERT_THAT(Li2.motifs_[2].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(Li2.motifs_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, LiH) {

    ASSERT_THAT(LiH.motifs_[0].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(LiH.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(LiH.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(LiH.motifs_[1].electrons_.indices(), ElementsAre());
    ASSERT_THAT(LiH.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(LiH.motifs_[1].type(), MotifType::Core);

    ASSERT_THAT(LiH.motifs_[2].electrons_.indices(), ElementsAre(2, 3));
    ASSERT_THAT(LiH.motifs_[2].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(LiH.motifs_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, H2covalent) {

    ASSERT_THAT(H2covalent.motifs_[0].electrons_.indices(), ElementsAre());
    ASSERT_THAT(H2covalent.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(H2covalent.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(H2covalent.motifs_[1].electrons_.indices(), ElementsAre());
    ASSERT_THAT(H2covalent.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(H2covalent.motifs_[1].type(), MotifType::Core);

    ASSERT_THAT(H2covalent.motifs_[2].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(H2covalent.motifs_[2].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(H2covalent.motifs_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, H2covalentMotifMerge) {

    auto H2motifs = H2covalent;

    ASSERT_EQ(H2motifs.motifs_.size(), 3);

    H2motifs.mergeMotifs({0, 1});

    ASSERT_EQ(H2motifs.motifs_.size(), 2);

    ASSERT_THAT(H2motifs.motifs_[0].electrons_.indices(), ElementsAre());
    ASSERT_THAT(H2motifs.motifs_[0].nuclei_.indices(), ElementsAre(0,1));
    ASSERT_EQ(H2motifs.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(H2motifs.motifs_[1].electrons_.indices(), ElementsAre(0,1));
    ASSERT_THAT(H2motifs.motifs_[1].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(H2motifs.motifs_[1].type(), MotifType::Valence);

    H2motifs.mergeMotifs({0, 1});
    ASSERT_EQ(H2motifs.motifs_.size(), 1);

    ASSERT_THAT(H2motifs.motifs_[0].electrons_.indices(), ElementsAre(0,1));
    ASSERT_THAT(H2motifs.motifs_[0].nuclei_.indices(), ElementsAre(0,1));
    ASSERT_EQ(H2motifs.motifs_[0].type(), MotifType::CoreValence);

}

TEST_F(AMotifsTest, H2ionicRight) {

    ASSERT_THAT(H2ionicRight.motifs_[0].electrons_.indices(), ElementsAre());
    ASSERT_THAT(H2ionicRight.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(H2ionicRight.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(H2ionicRight.motifs_[1].electrons_.indices(), ElementsAre());
    ASSERT_THAT(H2ionicRight.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(H2ionicRight.motifs_[1].type(), MotifType::Core);

    ASSERT_THAT(H2ionicRight.motifs_[2].electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(H2ionicRight.motifs_[2].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(H2ionicRight.motifs_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, BH3covalent) {
    ASSERT_THAT(BH3covalent.motifs_[0].electrons_.indices(), ElementsAre(0, 4));
    ASSERT_THAT(BH3covalent.motifs_[0].nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(BH3covalent.motifs_[0].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifs_[1].electrons_.indices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifs_[1].nuclei_.indices(), ElementsAre(1));
    ASSERT_EQ(BH3covalent.motifs_[1].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifs_[2].electrons_.indices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifs_[2].nuclei_.indices(), ElementsAre(2));
    ASSERT_EQ(BH3covalent.motifs_[2].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifs_[3].electrons_.indices(), ElementsAre());
    ASSERT_THAT(BH3covalent.motifs_[3].nuclei_.indices(), ElementsAre(3));
    ASSERT_EQ(BH3covalent.motifs_[3].type(), MotifType::Core);

    ASSERT_THAT(BH3covalent.motifs_[4].electrons_.indices(), ElementsAre(1, 5));
    ASSERT_THAT(BH3covalent.motifs_[4].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifs_[4].type(), MotifType::Valence);

    ASSERT_THAT(BH3covalent.motifs_[5].electrons_.indices(), ElementsAre(2, 6));
    ASSERT_THAT(BH3covalent.motifs_[5].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifs_[5].type(), MotifType::Valence);

    ASSERT_THAT(BH3covalent.motifs_[6].electrons_.indices(), ElementsAre(3, 7));
    ASSERT_THAT(BH3covalent.motifs_[6].nuclei_.indices(), ElementsAre());
    ASSERT_EQ(BH3covalent.motifs_[6].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, YAML) {
    auto node = YAML::convert<Motif>::encode(LiH.motifs_[0]);

    Motif decodedMotif;
    YAML::convert<Motif>::decode(node, decodedMotif);

    ASSERT_THAT(decodedMotif.electrons_.indices(), ElementsAre(0, 1));
    ASSERT_THAT(decodedMotif.nuclei_.indices(), ElementsAre(0));
    ASSERT_EQ(decodedMotif.type(), MotifType::Core);
}
