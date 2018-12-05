//
// Created by Michael Heuer on 08.05.18.
//

#include <gmock/gmock.h>
#include "MolecularGeometry.h"
#include "TestMolecules.h"

using namespace testing;

TEST(AMolecularGeometryTest, CorrectBraceInitialization) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.numberOfEntities(),4);
    ASSERT_EQ(molecule.atoms().numberOfEntities(),2);
    ASSERT_EQ(molecule.electrons().numberOfEntities(),2);
}

TEST(AMolecularGeometryTest, NumberedTypeByIndex) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.findNumberedTypeByIndex(0),NumberedElement(Element::H,0).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(1),NumberedElement(Element::H,1).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(2),NumberedSpin(Spin::alpha,0).toIntType());
    ASSERT_EQ(molecule.findNumberedTypeByIndex(3),NumberedSpin(Spin::beta,0).toIntType());
}

TEST(AMolecularGeometryTest, IndexFromNumberedType) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedElement(Element::H, 0).toIntType()).second, 0);
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedElement(Element::H, 1).toIntType()).second, 1);
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedSpin(Spin::alpha, 0).toIntType()).second, 2);//TODO ! IS THIS WANTED?
    ASSERT_EQ(molecule.findIndexByNumberedType(NumberedSpin(Spin::beta, 0).toIntType()).second, 3); // TODO ! IS THIS WANTED?
}