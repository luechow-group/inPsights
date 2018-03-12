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

        electronsVector1;
        electronsVector1.append(electron1);
        electronsVector1.append(electron2);
        electronsVector1.append(electron3);

        electronsVector2;
        electronsVector2.append(electron2);
        electronsVector2.append(electron1);
        electronsVector2.append(electron3);

        electronsVector3;
        electronsVector3.append(electron1);
        electronsVector3.append(electron3);
        electronsVector3.append(electron2);
    }
    ElectronsVector electronsVector1,electronsVector2,electronsVector3;
};

TEST_F(AElectronsVectorCollectionTest, InitialNumberEntities){

    ElectronsVectorCollection electronsVectorCollection;
    ASSERT_EQ(electronsVectorCollection.numberOfEntities(),0);
    ASSERT_EQ(electronsVectorCollection.positionsVectorCollection().numberOfEntities(),0);
    ASSERT_EQ(electronsVectorCollection.spinTypesVector().numberOfEntities(),0);
}

TEST_F(AElectronsVectorCollectionTest, NumberEntities){

    ElectronsVectorCollection electronsVectorCollection;

    electronsVectorCollection.append(electronsVector1);
    electronsVectorCollection.append(electronsVector1);

    ASSERT_EQ(electronsVectorCollection.numberOfEntities(),2);
    ASSERT_EQ(electronsVectorCollection[0].numberOfEntities(),3);
}


