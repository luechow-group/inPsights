//
// Created by Michael Heuer on 08.05.18.
//

#include <gmock/gmock.h>
#include "ParticleKit.h"
#include <TestMolecules.h>
#include <iomanip>

using namespace testing;

class AParticleKitTest : public ::testing::Test {
public:
    MolecularGeometry molecularGeometry;
    void SetUp() override {
        Particle<Spin > e1 = {Spin::alpha,{1, 2, 3}};
        Particle<Spin > e2 = {Spin::alpha,{1, 2, 3}};
        Particle<Spin > e3 = {Spin::beta ,{4, 5, 6}};
        molecularGeometry.electrons().append(e1);
        molecularGeometry.electrons().append(e2);
        molecularGeometry.electrons().append(e3);

        Particle<Element> a1 = {Element::He,{1, 2, 3}};
        Particle<Element> a2 = {Element::H ,{4, 5, 6}};
        molecularGeometry.atoms().append(a1);
        molecularGeometry.atoms().append(a2);
    };
};

TEST_F(AParticleKitTest, Constructor1) {
    ParticleKit::create(molecularGeometry.atoms(),0,2);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Element::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Element::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticleKitTest, Constructor2) {
    ParticleKit::create(molecularGeometry);

    // tests order
    ASSERT_EQ(ParticleKit::atomKit[0].first, Element::H);
    ASSERT_EQ(ParticleKit::atomKit[1].first, Element::He);

    ASSERT_EQ(ParticleKit::atomKit[0].second, 1);
    ASSERT_EQ(ParticleKit::atomKit[1].second, 1);

    ASSERT_EQ(ParticleKit::electronKit.first, 2);
    ASSERT_EQ(ParticleKit::electronKit.second, 1);
}

TEST_F(AParticleKitTest, WrongMultiplicityOrCharge) {
    EXPECT_DEATH(ParticleKit::create(molecularGeometry.atoms()),"");
}

TEST_F(AParticleKitTest, isSubsetQTrue) {

    AtomsVector atomsSubset;
    Atom a1 = {Element::He,{1, 2, 3}};
    atomsSubset.append(a1);

    ParticleKit::create(molecularGeometry.atoms(),0,2);

    ASSERT_TRUE(ParticleKit::isSubsetQ(atomsSubset));
}


TEST_F(AParticleKitTest, isSubsetQFalse) {
    MolecularGeometry mol(
            AtomsVector({{Element::H,{0,0,0}}}),
            ElectronsVector({{Spin::alpha,{0,0,0}}}));
    ParticleKit::create(mol);

    AtomsVector wrongElement({{Element::Ba,{1, 2, 3}}});
    ASSERT_FALSE(ParticleKit::isSubsetQ(wrongElement));

    AtomsVector wrongNumber({{Element::H,{0,0,0}},
                             {Element::H,{0,0,0}}});
    ASSERT_FALSE(ParticleKit::isSubsetQ(wrongNumber));

    MolecularGeometry mol2(
            AtomsVector({{Element::H,{0,0,0}}}),
            ElectronsVector({{Spin::beta,{0,0,0}}}));
    ASSERT_FALSE(ParticleKit::isSubsetQ(mol2));

    MolecularGeometry mol3(
            AtomsVector({{Element::H,{0,0,0}}}),
            ElectronsVector({{Spin::alpha,{0,0,0}},
                             {Spin::alpha,{0,0,0}}}));
    ASSERT_FALSE(ParticleKit::isSubsetQ(mol2));
}


TEST_F(AParticleKitTest, EnumeratedNumberedType) {
    ParticleKit::create({{Element::H,2},{Element::Ca,2},{Element::He,2}},{2,3});

    ASSERT_EQ(ParticleKit::getEnumeratedElementByIndex(0), EnumeratedElement (Element::H,0));
    ASSERT_EQ(ParticleKit::getEnumeratedElementByIndex(1), EnumeratedElement (Element::H,1));
    ASSERT_EQ(ParticleKit::getEnumeratedElementByIndex(4), EnumeratedElement (Element::He,0));
    ASSERT_EQ(ParticleKit::getEnumeratedElementByIndex(5), EnumeratedElement (Element::He,1));

    ASSERT_EQ(ParticleKit::getEnumeratedSpinByIndex(1), EnumeratedSpin(Spin::alpha,1));
    ASSERT_EQ(ParticleKit::getEnumeratedSpinByIndex(4), EnumeratedSpin(Spin::beta,2));

    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(0), EnumeratedType<int>(int(Element::H),0));
    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(1), EnumeratedType<int>(int(Element::H),1));
    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(4), EnumeratedType<int>(int(Element::He),0));
    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(5), EnumeratedType<int>(int(Element::He),1));
    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(7), EnumeratedType<int>(int(Spin::alpha),1));
    ASSERT_EQ(ParticleKit::getEnumeratedTypeByIndex(10), EnumeratedType<int>(int(Spin::beta),2));
}

TEST_F(AParticleKitTest, toString) {
    ParticleKit::create(TestMolecules::HeH::ElectronsInCores::normal);
    std::string ref = "ParticleKit:\n"
                      "------------\n"
                      "1*H, 1*He, 2*ea, 1*eb\n";
    ASSERT_EQ(ParticleKit::toString(),ref);
}

TEST_F(AParticleKitTest, electronsToKitPermutation) {
    ParticleKit::create(TestMolecules::threeElectrons::spinFlipped);
    auto original = TestMolecules::threeElectrons::spinFlipped.electrons();
    auto copy = original;

    //permute so that the result matches the kit
    copy.permute(ParticleKit::toKitPermutation(original));
    ASSERT_EQ(copy.typesVector(), ParticleKit::toSpinTypesVector());

    // permute copy back to the original order
    copy.permute(ParticleKit::fromKitPermutation(original));
    ASSERT_EQ(copy, original);
}

TEST_F(AParticleKitTest, CombinedPermutations) {

    ParticleKit::create(TestMolecules::threeElectrons::spinFlipped);
    auto original = TestMolecules::threeElectrons::spinFlipped.electrons();
    auto ref = original;
    auto test = original;

    Eigen::VectorXi indices(original.numberOfEntities());
    indices << 2,0,1;
    auto myperm = Eigen::PermutationMatrix<Eigen::Dynamic>(indices);
    auto tokit = ParticleKit::toKitPermutation(original);
    auto fromkit = ParticleKit::fromKitPermutation(original);

    // to kit
    ref.permute(tokit);
    test.permute(tokit);

    // both perms in sequence
    ref.permute(myperm);
    ref.permute(fromkit);

    // combined perm
    test.permute(fromkit*myperm);

    std::cout << std::endl;
    std::cout << original << std::endl;
    std::cout << ref << std::endl;
    std::cout << test << std::endl;
}

TEST_F(AParticleKitTest, atomsToKitPermutation) {
    ParticleKit::create(TestMolecules::CO2::nucleiPermuted);
    auto original = TestMolecules::CO2::nucleiPermuted.atoms();
    auto copy = original;

    //permute so that the result matches the kit
    copy.permute(ParticleKit::toKitPermutation(original));
    ASSERT_EQ(copy.typesVector(), ParticleKit::toElementTypesVector());

    // permute copy back to the original order
    copy.permute(ParticleKit::fromKitPermutation(original));
    ASSERT_EQ(copy, original);
}
