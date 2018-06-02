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

TEST_F(AMolecularGeometryTest, NumberedTypeByIndex) {
    ASSERT_EQ(molecule.findNumberedTypeByIndex(0),NumberedElement(Element::H,0).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(1),NumberedElement(Element::H,1).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(2),NumberedSpin(Spin::alpha,0).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(3),NumberedSpin(Spin::beta,0).toIntType());
}

TEST_F(AMolecularGeometryTest, IndexFromNumberedType) {
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedElement(Element::H, 0).toIntType()).second, 0);
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedElement(Element::H, 1).toIntType()).second, 1);
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedSpin(Spin::alpha, 0).toIntType()).second, 2);//TODO ! IS THIS WANTED?
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedSpin(Spin::beta, 0).toIntType()).second, 3); // TODO ! IS THIS WANTED?
}