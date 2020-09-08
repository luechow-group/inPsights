// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <ParticlesVector.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class AParticlesVectorTest : public Test {
public:

    Electron e0,e1,e2;
    Atom a0,a1;
    ElectronsVector electrons;
    AtomsVector atoms;

    void SetUp() override {
        e0 = {Spin::alpha,{1, 2, 3}};
        e1 = {Spin::alpha,{4, 5, 6}};
        e2 = {Spin::beta ,{7, 8, 9}};
        electrons.append(e0);
        electrons.append(e1);
        electrons.append(e2);

        a0 = {Element::H ,{1, 2, 3}};
        a1 = {Element::Og,{4, 5, 6}};
        atoms.append(a0);
        atoms.append(a1);
    };
};

TEST_F(AParticlesVectorTest, BraceInitialization) {

    ParticlesVector<Element> particlesVector(
            {{Element::H , {1, 2, 3}},
             {Element::Og, {4, 5, 6}}}
    );
}

TEST_F(AParticlesVectorTest, CopyConstructor) {

    auto electronsCopy = electrons;

    ASSERT_EQ(electronsCopy[0].type(),e0.type());
    ASSERT_EQ(electronsCopy[1].type(),e1.type());
    ASSERT_EQ(electronsCopy[2].type(),e2.type());

    ASSERT_EQ(electronsCopy[0].position(),e0.position());
    ASSERT_EQ(electronsCopy[1].position(),e1.position());
    ASSERT_EQ(electronsCopy[2].position(),e2.position());
}

TEST_F(AParticlesVectorTest, SpinTypeParticlesVector) {
    std::stringstream stringstream;
    stringstream << electrons;

    std::string expectedOutput = " 1 ea   1.00000   2.00000   3.00000\n"
                                 " 2 ea   4.00000   5.00000   6.00000\n"
                                 " 3 eb   7.00000   8.00000   9.00000\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

TEST_F(AParticlesVectorTest, ElementTypeParticlesVector) {
    std::stringstream stringstream;
    stringstream << atoms;

    std::string expectedOutput = " 1 H    1.00000   2.00000   3.00000\n"
                                 " 2 Og   4.00000   5.00000   6.00000\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

TEST_F(AParticlesVectorTest, CountTypeOccurence) {
    ASSERT_EQ(electrons.typesVector().countOccurence(Spin::alpha),2);
    ASSERT_EQ(atoms.typesVector().countOccurence(Element::H),1);
}

TEST_F(AParticlesVectorTest, Permute){
    auto e = electrons;

    VectorXi p(3);
    p << 2,0,1;

    e.permute(PermutationMatrix<Dynamic>(p));
    ASSERT_EQ(e[0].type(),e1.type());
    ASSERT_EQ(e[1].type(),e2.type());
    ASSERT_EQ(e[2].type(),e0.type());

    ASSERT_EQ(e[0].position(),e1.position());
    ASSERT_EQ(e[1].position(),e2.position());
    ASSERT_EQ(e[2].position(),e0.position());
}

/*
TEST_F(AParticlesVectorTest, IntegrationTest_PermuteAndTranslateAlphaElectrons){
    auto e = electrons;

    VectorXi p(2);
    p << 1,0;

    // swap positions
    e.positionsVector().slice({0,2},Reset::OnFinished).permute(PermutationMatrix<Dynamic>(p),Usage::NotFinished);
    e.positionsVector().translate({0,0,0.5},Usage::Finished);

    e.typesVector().slice({1,2}).permute(PermutationMatrix<Dynamic>(p));

    e.slice({0,2}).permute(PermutationMatrix<Dynamic>(p));

    ASSERT_EQ(e[0].type(),e2.type());
    ASSERT_EQ(e[1].type(),e0.type());
    ASSERT_EQ(e[2].type(),e1.type());

    ASSERT_EQ(e[0].position(),e0.position()+Vector3d(0,0,0.5));
    ASSERT_EQ(e[1].position(),e1.position()+Vector3d(0,0,0.5));
    ASSERT_EQ(e[2].position(),e2.position());
}*/

TEST_F(AParticlesVectorTest, LinkedParticles){
    ElectronsVector e = electrons;

    auto l0 = e.linkedParticle(0);
    auto l1 = e.linkedParticle(1);
    auto l2 = e.linkedParticle(2);

    ASSERT_EQ(e[0].type(),l0->type());
    ASSERT_EQ(e[1].type(),l1->type());
    ASSERT_EQ(e[2].type(),l2->type());

    ASSERT_EQ(e[0].position(),l0->position());
    ASSERT_EQ(e[1].position(),l1->position());
    ASSERT_EQ(e[2].position(),l2->position());

    auto changed = Eigen::Vector3d({1,1,1});
    l1->setPosition(changed);
    l2->setType(Spin::none);


    ASSERT_EQ(e[0].type(),l0->type());
    ASSERT_EQ(e[1].type(),l1->type());
    ASSERT_EQ(e[2].type(),Spin::none);

    ASSERT_EQ(e[0].position(),l0->position());
    ASSERT_EQ(e[1].position(),changed);
    ASSERT_EQ(e[2].position(),l2->position());
}

TEST_F(AParticlesVectorTest, AccessWithIndexList){
    std::list<long> indices({0,2});

    auto newVector = electrons[indices];
    ASSERT_EQ(newVector.numberOfEntities(), 2);
    ASSERT_EQ(newVector[0].position(), electrons[0].position());
    ASSERT_EQ(newVector[0].type(), electrons[0].type());
    ASSERT_EQ(newVector[1].position(), electrons[2].position());
    ASSERT_EQ(newVector[1].type(), electrons[2].type());
}

TEST_F(AParticlesVectorTest, AccessFirstElements){
    auto newVector = electrons.head(2);
    ASSERT_EQ(newVector.numberOfEntities(), 2);
    ASSERT_EQ(newVector[0].position(), electrons[0].position());
    ASSERT_EQ(newVector[0].type(), electrons[0].type());
    ASSERT_EQ(newVector[1].position(), electrons[1].position());
    ASSERT_EQ(newVector[1].type(), electrons[1].type());
}

TEST_F(AParticlesVectorTest, AccessLastElements){
    auto newVector = electrons.tail(2);
    ASSERT_EQ(newVector.numberOfEntities(), 2);
    ASSERT_EQ(newVector[0].position(), electrons[1].position());
    ASSERT_EQ(newVector[0].type(), electrons[1].type());
    ASSERT_EQ(newVector[1].position(), electrons[2].position());
    ASSERT_EQ(newVector[1].type(), electrons[2].type());
}
