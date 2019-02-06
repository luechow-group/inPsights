//
// Created by heuer on 06.02.19.
//

#include <gmock/gmock.h>
#include <MotifAnalysis.h>

using namespace testing;

class AMotifAnalysisTest : public ::testing::Test {
public:
    Eigen::MatrixXb A,B;
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
        // Two pairs: (0-1,) (2-3)
        A << 0, 1, 0, 0, \
             1, 0, 0, 0, \
             0, 0, 0, 1, \
             0, 0, 1, 0;
        B = Eigen::MatrixXb(4, 4);
        B << 0, 1, 1, 0, \
             1, 0, 1, 0, \
             1, 1, 0, 0, \
             0, 0, 0, 0;
    };
};

TEST_F(AMotifAnalysisTest, Test){

    MotifAnalysis::Motifs motifs(B, molecule);

    ASSERT_THAT(motifs.motifVector[0].electronIndices(), ElementsAre(0,1));
    ASSERT_EQ(motifs.motifVector[0].type(), MotifAnalysis::MotifType::Core);

    ASSERT_THAT(motifs.motifVector[1].electronIndices(), ElementsAre(2));
    ASSERT_EQ(motifs.motifVector[1].type(), MotifAnalysis::MotifType::Core);

    ASSERT_THAT(motifs.motifVector[2].electronIndices(), ElementsAre(3));
    ASSERT_EQ(motifs.motifVector[2].type(), MotifAnalysis::MotifType::Valence);

    //YAML::Emitter out;
    //out << motifs.motifVector;
    //std::cout << out.c_str() << std::endl;
}