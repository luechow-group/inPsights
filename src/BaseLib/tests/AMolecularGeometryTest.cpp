//
// Created by Michael Heuer on 08.05.18.
//

#include <gtest/gtest.h>
#include "MolecularGeometry.h"

using namespace testing;

class AMolecularGeometryTest : public Test {
public:
    MolecularGeometry molecule;
    void SetUp() override {
        molecule = {AtomsVector(
                            {{Elements::ElementType::H,{0,0, 0.3705}},
                             {Elements::ElementType::H,{0,0,-0.3705}}}),
                    ElectronsVector(
                            {{Spins::SpinType::alpha,{0,0, 0.3705}},
                             {Spins::SpinType::beta, {0,0,-0.3705}}})};
    }

};

TEST_F(AMolecularGeometryTest, CorrectBraceInitialization) {
    ASSERT_EQ(molecule.numberOfEntities(),4);
    ASSERT_EQ(molecule.atoms().numberOfEntities(),2);
    ASSERT_EQ(molecule.electrons().numberOfEntities(),2);
}
