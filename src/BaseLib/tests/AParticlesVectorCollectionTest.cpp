// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <ParticlesVectorCollection.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class AParticlesVectorCollectionTest : public Test {
public:

    Eigen::Vector3d pos1{1,2,3};
    Eigen::Vector3d pos2{4,5,6};
    ElectronsVector electronsVector;
    ElectronsVectorCollection electronsVectorCollection;
    AtomsVector atomsVector;
    AtomsVectorCollection atomsVectorCollection;

    void SetUp() override {
        Particle<Spin > e1 = {Spin::alpha,pos1};
        Particle<Spin > e2 = {Spin::beta,pos2};
        electronsVector.append(e1);
        electronsVector.append(e2);

        electronsVectorCollection.append(electronsVector);
        electronsVectorCollection.append(electronsVector);

        Particle<Element> a1 = {Element::H ,pos1};
        Particle<Element> a2 = {Element::Og,pos2};
        atomsVector.append(a1);
        atomsVector.append(a2);
        atomsVectorCollection.append(atomsVector);
        atomsVectorCollection.append(atomsVector);
    };
};

/*TEST_F(AParticlesVectorCollectionTest, Constructor) {
    EXPECT_TRUE(false);
}

TEST_F(AParticlesVectorCollectionTest, CopyConstructor) {
    EXPECT_TRUE(false);
}*/

TEST_F(AParticlesVectorCollectionTest, SpinTypeParticlesVectorCollection) {
    std::stringstream stringstream;
    stringstream << electronsVectorCollection;

    std::string expectedOutput = "Vector 1:\n"
                                 " 1 ea   1.00000   2.00000   3.00000\n"
                                 " 2 eb   4.00000   5.00000   6.00000\n"
                                 "\n"
                                 "Vector 2:\n"
                                 " 1 ea   1.00000   2.00000   3.00000\n"
                                 " 2 eb   4.00000   5.00000   6.00000\n"
                                 "\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

TEST_F(AParticlesVectorCollectionTest, ElementTypeParticlesVectorCollection) {
    std::stringstream stringstream;
    stringstream << atomsVectorCollection;

    std::string expectedOutput = "Vector 1:\n"
                                 " 1 H    1.00000   2.00000   3.00000\n"
                                 " 2 Og   4.00000   5.00000   6.00000\n"
                                 "\n"
                                 "Vector 2:\n"
                                 " 1 H    1.00000   2.00000   3.00000\n"
                                 " 2 Og   4.00000   5.00000   6.00000\n"
                                 "\n";
    ASSERT_EQ(stringstream.str(), expectedOutput);
}

