//
// Created by Michael Heuer on 08.05.18.
//

#include <gmock/gmock.h>
#include "MolecularGeometry.h"
#include "TestMolecules.h"

using namespace testing;

class AMolecularGeometryTest: public Test {
public:
    MolecularGeometry mol;
    void SetUp() override {

        mol = {
                AtomsVector({
                    {Element::H ,{ 1, 0, 0}},
                    {Element::He,{-1, 0, 0}}
                }),
                ElectronsVector({
                    {Spin::alpha,{   1, 0, 0}},
                    {Spin::alpha,{   1, 0, 0}},
                    {Spin::alpha,{ 1.1, 0, 0}},
                    {Spin::beta,{   -1, 0, 0}},
                    {Spin::beta,{-1.09, 0, 0}}
                })};
    };
};

TEST_F(AMolecularGeometryTest, CorrectBraceInitialization) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.numberOfEntities(),4);
    ASSERT_EQ(molecule.atoms().numberOfEntities(),2);
    ASSERT_EQ(molecule.electrons().numberOfEntities(),2);
}

TEST_F(AMolecularGeometryTest, EnumeratedTypeByIndex) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(0),NumberedElement(Element::H,0).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(1),NumberedElement(Element::H,1).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(2),NumberedSpin(Spin::alpha,0).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(3),NumberedSpin(Spin::beta,0).toIntType());
}

TEST_F(AMolecularGeometryTest, IndexFromEnumeratedType) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.findIndexByEnumeratedType(NumberedElement(Element::H, 0).toIntType()).second, 0);
    ASSERT_EQ(molecule.findIndexByEnumeratedType(NumberedElement(Element::H, 1).toIntType()).second, 1);
    ASSERT_EQ(molecule.findIndexByEnumeratedType(NumberedSpin(Spin::alpha, 0).toIntType()).second, 2);//TODO ! IS THIS WANTED?
    ASSERT_EQ(molecule.findIndexByEnumeratedType(NumberedSpin(Spin::beta, 0).toIntType()).second, 3); // TODO ! IS THIS WANTED?
}

TEST_F(AMolecularGeometryTest, CoreElectrons) {
    double thresh = 0.1;
    ASSERT_THAT(mol.coreElectronsIndices(0,thresh), ElementsAre(0,1));
    ASSERT_THAT(mol.coreElectronsIndices(1,thresh), ElementsAre(3,4));
    ASSERT_THAT(mol.coreElectronsIndices(thresh), ElementsAre(0,1,3,4));
}

TEST_F(AMolecularGeometryTest, ValenceElectrons) {
    double thresh = 0.1;
    ASSERT_THAT(mol.valenceElectronsIndices(thresh), ElementsAre(2));
}
