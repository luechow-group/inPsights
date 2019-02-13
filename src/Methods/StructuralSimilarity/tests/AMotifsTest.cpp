//
// Created by heuer on 06.02.19.
//

#include <gmock/gmock.h>
#include <Motifs.h>

using namespace testing;

class AMotifsTest : public ::testing::Test {
public:
    Eigen::MatrixXb A;
    MolecularGeometry molecule;
    void SetUp() override {

        molecule = {
                AtomsVector({{Element::H, {0, 0, 1}},
                             {Element::H, {0, 0,-1}}}),
                ElectronsVector({{Spin::alpha, {0, 0, 1}},
                                 {Spin::beta,  {0, 0, 1}},
                                 {Spin::alpha, {0, 0,-1}},
                                 {Spin::beta,  {0, 0,-1.1}},
                })};

        A = Eigen::MatrixXb(4, 4);
        A << 0, 1, 1, 0, \
             1, 0, 1, 0, \
             1, 1, 0, 0, \
             0, 0, 0, 0;
    };
};

TEST_F(AMotifsTest, Assignment) {
    Motifs motifs(A, molecule);

    ASSERT_THAT(motifs.motifVector_[0].electronIndices(), ElementsAre(0,1));
    ASSERT_THAT(motifs.motifVector_[0].atomIndices(), ElementsAre(0));
    ASSERT_EQ(motifs.motifVector_[0].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[1].electronIndices(), ElementsAre(2));
    ASSERT_THAT(motifs.motifVector_[1].atomIndices(), ElementsAre(1));
    ASSERT_EQ(motifs.motifVector_[1].type(), MotifType::Core);

    ASSERT_THAT(motifs.motifVector_[2].electronIndices(), ElementsAre(3));
    ASSERT_THAT(motifs.motifVector_[2].atomIndices(), ElementsAre());
    ASSERT_EQ(motifs.motifVector_[2].type(), MotifType::Valence);
}
