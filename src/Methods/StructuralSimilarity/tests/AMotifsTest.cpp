//
// Created by heuer on 06.02.19.
//

#include <gmock/gmock.h>
#include <Motifs.h>

using namespace testing;

class AMotifsTest : public ::testing::Test {
public:
    Eigen::MatrixXb A;

    Motifs motifs;
    void SetUp() override {

        MolecularGeometry molecule = {
                AtomsVector({{Element::Li, {0, 0, 1}},
                             {Element::H, {0, 0,-1}}}),
                ElectronsVector({{Spin::alpha, {0, 0, 1}},
                                 {Spin::beta,  {0, 0, 1}},
                                 {Spin::alpha, {0, 0,-1}},
                                 {Spin::beta,  {0, 0,-1.1}},
                })};


        A = Eigen::MatrixXb(4, 4);
        A << 0, 1, 0, 0, \
             1, 0, 0, 0, \
             0, 0, 0, 1, \
             0, 0, 1, 0;

        motifs = Motifs(A, molecule);
    };
};

TEST_F(AMotifsTest, Sorting) {

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre());
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[2].electronIndices(), ElementsAre(2,3));
    ASSERT_THAT(motifs.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(motifs.motifVector_[2].type(), MotifType::Valence);
}

TEST_F(AMotifsTest, YAML) {
    auto node = YAML::convert<Motif>::encode(motifs.motifVector_[0]);

    Motif decodedMotif;
    YAML::convert<Motif>::decode(node, decodedMotif);

    ASSERT_THAT(decodedMotif.electronIndices(),ElementsAre(0,1));
    ASSERT_THAT(decodedMotif.atomIndices(),ElementsAre(0));
    ASSERT_EQ(decodedMotif.type(), MotifType::Core);
}
