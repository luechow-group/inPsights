//
// Created by Leonard Reuter on 09.03.18.
//

#include <gtest/gtest.h>
#include "AtomsVectorCollection.h"

using namespace testing;
using namespace Eigen;

class AAtomsVectorCollectionTest : public Test {
public:
    void SetUp() override {
        Atom atom1(Vector3d(1,2,3),Elements::ElementType::Ag);
        Atom atom2(Vector3d(-1,0,3.5),Elements::ElementType::Au);
        Atom atom3(Vector3d(-7.3,0.5,9),Elements::ElementType::C);

        atomCollection1;
        atomCollection1.append(atom1);
        atomCollection1.append(atom2);
        atomCollection1.append(atom3);

        atomCollection2;
        atomCollection2.append(atom2);
        atomCollection2.append(atom1);
        atomCollection2.append(atom3);

        atomCollection3;
        atomCollection3.append(atom1);
        atomCollection3.append(atom3);
        atomCollection3.append(atom2);
    }
    AtomCollection atomCollection1,atomCollection2,atomCollection3;
};

TEST_F(AAtomsVectorCollectionTest, InitialNumberEntities){

    AtomsVectorCollection atomsVectorCollection;
    ASSERT_EQ(atomsVectorCollection.numberOfEntities(),0);
    ASSERT_EQ(atomsVectorCollection.positionsVectorCollection().numberOfEntities(),0);
    ASSERT_EQ(atomsVectorCollection.elementTypeCollection().numberOfEntities(),0);
}

TEST_F(AAtomsVectorCollectionTest, NumberEntities){

    AtomsVectorCollection atomsVectorCollection;

    atomsVectorCollection.append(atomCollection1);
    atomsVectorCollection.append(atomCollection1);

    ASSERT_EQ(atomsVectorCollection.numberOfEntities(),2);
    ASSERT_EQ(atomsVectorCollection[0].numberOfEntities(),3);
}


