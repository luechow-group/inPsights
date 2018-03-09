//
// Created by Leonard Reuter on 09.03.18.
//

#include <gtest/gtest.h>
#include "ElectronsVectorCollection.h"

using namespace testing;
using namespace Eigen;

class AElectronsVectorCollectionTest : public Test {
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

TEST_F(AElectronsVectorCollectionTest, InitialNumberEntities){

    ElectronsVectorCollection electronsVectorCollection;
    ASSERT_EQ(electronsVectorCollection.numberOfEntities(),0);
    ASSERT_EQ(electronsVectorCollection.positionsVectorCollection().numberOfEntities(),0);
    ASSERT_EQ(electronsVectorCollection.spinTypesVector().numberOfEntities(),0);
}

TEST_F(AElectronsVectorCollectionTest, NumberEntities){

    ElectronsVectorCollection electronsVectorCollection;

    electronsVectorCollection.append(electronCollection1);
    electronsVectorCollection.append(electronCollection1);

    ASSERT_EQ(electronsVectorCollection.numberOfEntities(),2);
    ASSERT_EQ(electronsVectorCollection[0].numberOfEntities(),3);
}


