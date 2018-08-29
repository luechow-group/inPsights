//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "PositionsVector.h"

using namespace testing;
using namespace Eigen;

class APositionsVectorTest : public Test {
public:
    void SetUp() override {
        position0 = Vector3d(1,2,3);
        position1 = Vector3d(4,5,6);
        position2 = Vector3d(7,8,9);
    }
    Vector3d position0,position1,position2;
};

TEST_F(APositionsVectorTest, NumberOfEntities){

    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    ASSERT_EQ(positionsVector.numberOfEntities(),3);
}

TEST_F(APositionsVectorTest, IndexOperator){

    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
    ASSERT_EQ(positionsVector[-1],position2);
    ASSERT_EQ(positionsVector[-2],position1);
    ASSERT_EQ(positionsVector[-3],position0);

    EXPECT_DEATH(positionsVector[3],"");
    EXPECT_DEATH(positionsVector[-4],"");
}

TEST_F(APositionsVectorTest, insert){
    VectorXd positions(6);
    positions<< position0, position2;
    Vector3d position(position1);
    PositionsVector positionsVector(positions);

    positionsVector.insert(position,1);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}


TEST_F(APositionsVectorTest, prepend) {
    VectorXd positions(6);
    positions<< position1, position2;
    Vector3d position(position0);
    PositionsVector positionsVector(positions);

    positionsVector.prepend(position);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}

TEST_F(APositionsVectorTest, append) {
    VectorXd positions(6);
    positions<< position0, position1;
    Vector3d position(position2);
    PositionsVector positionsVector(positions);

    positionsVector.append(position);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}

TEST_F(APositionsVectorTest, print) {
    VectorXd positions(9);
    Eigen::Vector3d a(0,-2.01238,7.13);
    Eigen::Vector3d b(1e-10,-9.8987,-12.238);
    Eigen::Vector3d c(-99.9,99.9,9.99);
    positions<< a, b, c;
    PositionsVector positionsVector(positions);
    std::ostringstream os;
    os << positionsVector;
    ASSERT_EQ(os.str()," 1    0.00000  -2.01238   7.13000\n 2    0.00000  -9.89870 -12.23800\n 3  -99.90000  99.90000   9.99000\n");
}

TEST_F(APositionsVectorTest, permute) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    Eigen::VectorXi p(3);
    p << 0,2,1;
    PermutationMatrix<Eigen::Dynamic> perm(p);
    positionsVector.permute(perm);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position2);
    ASSERT_EQ(positionsVector[2],position1);
}

TEST_F(APositionsVectorTest, slice) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    ASSERT_EQ(*positionsVector.entity(0).positionsRefPtr_, positions.segment(0,3));
    ASSERT_EQ(*positionsVector.entity(1).positionsRefPtr_, positions.segment(3,3));
    ASSERT_EQ(*positionsVector.entity(2).positionsRefPtr_, positions.segment(6,3));

    ASSERT_EQ(*positionsVector.slice({0,1}).positionsRefPtr_,*positionsVector.entity(0).positionsRefPtr_);
    ASSERT_EQ(*positionsVector.slice({1,1}).positionsRefPtr_,*positionsVector.entity(1).positionsRefPtr_);
    ASSERT_EQ(*positionsVector.slice({2,1}).positionsRefPtr_,*positionsVector.entity(2).positionsRefPtr_);

    ASSERT_EQ(*positionsVector.slice({0,2}).positionsRefPtr_, positions.segment(0,6));
    ASSERT_EQ(*positionsVector.slice({1,2}).positionsRefPtr_, positions.segment(3,6));
    ASSERT_EQ(*positionsVector.slice({0,3}).positionsRefPtr_, positions.segment(0,9));

    ASSERT_EQ(*positionsVector.slice({0,3}).positionsRefPtr_,*positionsVector.all().positionsRefPtr_);

    EXPECT_DEATH(*positionsVector.entity(-1).positionsRefPtr_,"");
    EXPECT_DEATH(*positionsVector.slice({0,4}).positionsRefPtr_,"");
}

TEST_F(APositionsVectorTest, translate) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    VectorXd expected[2] = {VectorXd(9),VectorXd(9)};
    expected[0] << 1,2,3,5,6,7,7,8,9;
    expected[1] << 1,2,3,6,7,8,8,9,10;

    positionsVector.entity(1).translate({1,1,1});
    std::cout << positionsVector << std::endl;
    ASSERT_EQ(positionsVector.positionsAsEigenVector(),expected[0]);
    positionsVector.slice({1,2}).translate({1,1,1});
    std::cout << positionsVector << std::endl;
    ASSERT_EQ(positionsVector.positionsAsEigenVector(),expected[1]);
}