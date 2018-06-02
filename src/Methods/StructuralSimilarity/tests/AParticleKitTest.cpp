//
// Created by Michael Heuer on 08.05.18.
//

#include <gtest/gtest.h>
#include "ParticleKit.h"
#include <iomanip>

using namespace testing;

class AParticleKitTest : public ::testing::Test {
public:
    MolecularGeometry molecularGeometry;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e2 = {Spins::SpinType::alpha,{1, 2, 3}};
        Particle<Spins::SpinType > e3 = {Spins::SpinType::beta ,{4, 5, 6}};
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
}


TEST_F(AParticleKitTest, NumberedType) {
    ParticleKit::create({{Element::H,2},{Element::Ca,2},{Element::He,2}},{2,3});

    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(0), NumberedElement (Element::H,0));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(1), NumberedElement (Element::H,1));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(4), NumberedElement (Element::He,0));
    ASSERT_EQ(ParticleKit::getNumberedElementByIndex(5), NumberedElement (Element::He,1));

    ASSERT_EQ(ParticleKit::getNumberedSpinByIndex(1), NumberedSpin(Spin::alpha,1));
    ASSERT_EQ(ParticleKit::getNumberedSpinByIndex(4), NumberedSpin(Spin::beta,2));



    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(0), NumberedType<int>(int(Element::H),0));
    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(1), NumberedType<int>(int(Element::H),1));
    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(4), NumberedType<int>(int(Element::He),0));
    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(5), NumberedType<int>(int(Element::He),1));
    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(7), NumberedType<int>(int(Spin::alpha),1));
    ASSERT_EQ(ParticleKit::getNumberedTypeByIndex(10), NumberedType<int>(int(Spin::beta),2));
}
