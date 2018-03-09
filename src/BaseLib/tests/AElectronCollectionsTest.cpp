//
// Created by Leonard Reuter on 09.03.18.
//

#include <gtest/gtest.h>
#include "ElectronCollections.h"

using namespace testing;
using namespace Eigen;

class AElectronCollectionsTest : public Test {
public:
    void SetUp() override {
        Electron electron1(Vector3d(1,2,3),Spin::SpinType::alpha);
        Electron electron2(Vector3d(-1,0,3.5),Spin::SpinType::beta);
        Electron electron3(Vector3d(-7.3,0.5,9),Spin::SpinType::none);

        electronCollection1;
        electronCollection1.append(electron1);
        electronCollection1.append(electron2);
        electronCollection1.append(electron3);

        electronCollection2;
        electronCollection2.append(electron2);
        electronCollection2.append(electron1);
        electronCollection2.append(electron3);

        electronCollection3;
        electronCollection3.append(electron1);
        electronCollection3.append(electron3);
        electronCollection3.append(electron2);
    }
    ElectronCollection electronCollection1,electronCollection2,electronCollection3;
};

TEST_F(AElectronCollectionsTest, InitialNumberEntities){

    ElectronCollections electronCollections;
    ASSERT_EQ(electronCollections.numberOfEntities(),0);
    ASSERT_EQ(electronCollections.positionCollections().numberOfEntities(),0);
    ASSERT_EQ(electronCollections.spinTypeCollection().numberOfEntities(),0);
}

TEST_F(AElectronCollectionsTest, NumberEntities){

    ElectronCollections electronCollections;

    electronCollections.append(electronCollection1);
    electronCollections.append(electronCollection1);

    ASSERT_EQ(electronCollections.numberOfEntities(),2);
    ASSERT_EQ(electronCollections[0].numberOfEntities(),3);
}


