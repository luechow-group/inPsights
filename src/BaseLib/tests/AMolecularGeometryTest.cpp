// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(0),EnumeratedElement(Element::H,0).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(1),EnumeratedElement(Element::H,1).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(2),EnumeratedSpin(Spin::alpha,0).toIntType());
    ASSERT_EQ(molecule.findEnumeratedTypeByIndex(3),EnumeratedSpin(Spin::beta,0).toIntType());
}

TEST_F(AMolecularGeometryTest, IndexFromEnumeratedType) {
    auto molecule = TestMolecules::H2::ElectronsInCores::normal;

    ASSERT_EQ(molecule.findIndexByEnumeratedType(EnumeratedElement(Element::H, 0).toIntType()).second, 0);
    ASSERT_EQ(molecule.findIndexByEnumeratedType(EnumeratedElement(Element::H, 1).toIntType()).second, 1);
    ASSERT_EQ(molecule.findIndexByEnumeratedType(EnumeratedSpin(Spin::alpha, 0).toIntType()).second, 2);//TODO ! IS THIS WANTED?
    ASSERT_EQ(molecule.findIndexByEnumeratedType(EnumeratedSpin(Spin::beta, 0).toIntType()).second, 3); // TODO ! IS THIS WANTED?
}

TEST_F(AMolecularGeometryTest, CoreElectrons) {
    double thresh = 0.1;
    ASSERT_THAT(mol.coreElectronsIndices(0,thresh), ElementsAre(0,1));
    ASSERT_THAT(mol.coreElectronsIndices(1,thresh), ElementsAre(3,4));
    ASSERT_THAT(mol.coreElectronsIndices(thresh), ElementsAre(0,1,3,4));
}

TEST_F(AMolecularGeometryTest, ValenceElectrons) {
    double thresh = 0.1;
    ASSERT_THAT(mol.nonCoreElectronsIndices(thresh), ElementsAre(2));
}

TEST_F(AMolecularGeometryTest, Positions) {
    auto mol = TestMolecules::H2::ElectronsInCores::normal;

    auto d = mol.atoms().positionsVector().entityLength();
    auto M = mol.atoms().numberOfEntities();
    auto N = mol.electrons().numberOfEntities();

    ASSERT_TRUE(
            mol.positions().asEigenVector().head(M*d).isApprox(mol.atoms().positionsVector().asEigenVector()));
    ASSERT_TRUE(
            mol.positions().asEigenVector().tail(N*d).isApprox(mol.electrons().positionsVector().asEigenVector()));
}

TEST_F(AMolecularGeometryTest, EqualOperator) {
    auto mol = TestMolecules::H2::ElectronsInCores::normal;
    auto sameMol = mol;
    auto otherMol = TestMolecules::H2::ElectronsInCores::flippedSpins;

    ASSERT_TRUE(mol == sameMol);
    ASSERT_FALSE(mol == otherMol);
    ASSERT_TRUE(mol != otherMol);
}

TEST_F(AMolecularGeometryTest, SplitAllParticlePermutation) {
    auto mol = TestMolecules::BH3::ionicMinimalRotatedPermuted;
    Eigen::VectorXi allParticlePerm(mol.numberOfEntities());
    allParticlePerm << 3, 0, 1, 2, 4, 6, 5, 7;

    auto molPerms = mol.splitAllParticlePermutation(Eigen::PermutationMatrix<Eigen::Dynamic>(allParticlePerm));
    Eigen::VectorXi expectedNuclearPermutation(mol.atoms().numberOfEntities());
    Eigen::VectorXi expectedElectronicPermutation(mol.electrons().numberOfEntities());
    expectedNuclearPermutation << 3,0,1,2;
    expectedElectronicPermutation << 0,2,1,3;

    ASSERT_EQ(molPerms.nuclearPermutation.indices(), expectedNuclearPermutation);
    ASSERT_EQ(molPerms.electronicPermutation.indices(), expectedElectronicPermutation);

}
