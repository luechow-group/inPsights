//
// Created by Michael Heuer on 23.05.18.
//
#include <gtest/gtest.h>
#include "StructuralSimilarity.h"

class AStructuralSimilarityTest : public ::testing::Test {
public:
    MolecularGeometry A,B;
    void SetUp() override {
        A = {AtomsVector(
                {{Elements::ElementType::H,{0,0, 0.37}},
                 {Elements::ElementType::H,{0,0,-0.37}}}),
             ElectronsVector(
                     {{Spins::SpinType::alpha,{0,0, 0.37}},
                      {Spins::SpinType::alpha,{0,0,-0.37}}})
        };

        B = {AtomsVector(
                {{Elements::ElementType::H,{0,0, 0.37}},
                 {Elements::ElementType::H,{0,0,-0.37}}}),
             ElectronsVector(
                     {{Spins::SpinType::alpha,{0,0, 0.37}},
                      {Spins::SpinType::alpha,{0,0,-0.37}}})
        };
    }
};

TEST_F(AStructuralSimilarityTest , Test) {
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;
    ParticleKit::create({{Elements::ElementType::H,2},{Elements::ElementType::He,2}},{2,2});

    ASSERT_TRUE(ParticleKit::isSubsetQ(A));
    ASSERT_TRUE(ParticleKit::isSubsetQ(B));

    double result = StructuralSimilarity::stucturalSimilarity(A,B,1);
    std::cout << "FINAL RESULT:" << result << std::endl;
    ASSERT_EQ(result, 1.0);
}