//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "PositionCollection.h"

using namespace testing;
using namespace Eigen;

class APositionCollectionTest : public Test {
public:
    void SetUp() override {
        position0 = Vector3d(1,2,3);
        position1 = Vector3d(4,5,6);
        position2 = Vector3d(7,8,9);
    }
    Vector3d position0,position1,position2;
};

TEST_F(APositionCollectionTest, NumberOfEntities){

    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionCollection positionCollection(positions);

    ASSERT_EQ(positionCollection.numberOfEntities(),3);
}

TEST_F(APositionCollectionTest, IndexOperator){

    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionCollection positionCollection(positions);

    ASSERT_EQ(positionCollection[0],position0);
    ASSERT_EQ(positionCollection[1],position1);
    ASSERT_EQ(positionCollection[2],position2);
    ASSERT_EQ(positionCollection[-1],position2);
    ASSERT_EQ(positionCollection[-2],position1);
    ASSERT_EQ(positionCollection[-3],position0);

    EXPECT_DEATH(positionCollection[3],"");
    EXPECT_DEATH(positionCollection[-4],"");
}

TEST_F(APositionCollectionTest, insert){
    VectorXd positions(6);
    positions<< position0, position2;
    Vector3d position(position1);
    PositionCollection positionCollection(positions);

    positionCollection.insert(position,1);

    ASSERT_EQ(positionCollection[0],position0);
    ASSERT_EQ(positionCollection[1],position1);
    ASSERT_EQ(positionCollection[2],position2);
}


TEST_F(APositionCollectionTest, prepend) {
    VectorXd positions(6);
    positions<< position1, position2;
    Vector3d position(position0);
    PositionCollection positionCollection(positions);

    positionCollection.prepend(position);

    ASSERT_EQ(positionCollection[0],position0);
    ASSERT_EQ(positionCollection[1],position1);
    ASSERT_EQ(positionCollection[2],position2);
}

TEST_F(APositionCollectionTest, append) {
    VectorXd positions(6);
    positions<< position0, position1;
    Vector3d position(position2);
    PositionCollection positionCollection(positions);

    positionCollection.append(position);

    ASSERT_EQ(positionCollection[0],position0);
    ASSERT_EQ(positionCollection[1],position1);
    ASSERT_EQ(positionCollection[2],position2);
}

TEST_F(APositionCollectionTest, print) {
    VectorXd positions(9);
    Eigen::Vector3d a(0,-2.01238,7.13);
    Eigen::Vector3d b(1e-10,-9.8987,-12.238);
    Eigen::Vector3d c(-99.9,99.9,9.99);
    positions<< a, b, c;
    PositionCollection positionCollection(positions);
    std::ostringstream os;
    os << positionCollection;
    ASSERT_EQ(os.str()," 1    0.00000  -2.01238   7.13000\n 2    0.00000  -9.89870 -12.23800\n 3  -99.90000  99.90000   9.99000\n");
}