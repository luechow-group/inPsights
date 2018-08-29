//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include <PositionsVector.h>
#include <NaturalConstants.h>

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

TEST_F(APositionsVectorTest, Insert){
    VectorXd positions(6);
    positions<< position0, position2;
    Vector3d position(position1);
    PositionsVector positionsVector(positions);

    positionsVector.insert(position,1);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}


TEST_F(APositionsVectorTest, Prepend) {
    VectorXd positions(6);
    positions<< position1, position2;
    Vector3d position(position0);
    PositionsVector positionsVector(positions);

    positionsVector.prepend(position);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}

TEST_F(APositionsVectorTest, Append) {
    VectorXd positions(6);
    positions<< position0, position1;
    Vector3d position(position2);
    PositionsVector positionsVector(positions);

    positionsVector.append(position);

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position2);
}

TEST_F(APositionsVectorTest, Print) {
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

TEST_F(APositionsVectorTest, PermuteAll) {
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

TEST_F(APositionsVectorTest, PermuteAllCyclic) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    Eigen::VectorXi p(3);
    p << 2,0,1;
    PermutationMatrix<Eigen::Dynamic> perm(p);
    positionsVector.permute(perm);

    ASSERT_EQ(positionsVector[0],position1);
    ASSERT_EQ(positionsVector[1],position2);
    ASSERT_EQ(positionsVector[2],position0);
}

TEST_F(APositionsVectorTest, PermuteSlice) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    Eigen::VectorXi p1(2);
    p1 << 1,0;
    positionsVector.slice({1,2}).permute(PermutationMatrix<Eigen::Dynamic>(p1));

    ASSERT_EQ(positionsVector[0],position0);
    ASSERT_EQ(positionsVector[1],position2);
    ASSERT_EQ(positionsVector[2],position1);

    Eigen::VectorXi p2(3);
    p2 << 2,0,1;

    // test resetToAll()
    positionsVector.permute(PermutationMatrix<Eigen::Dynamic>(p2));
    ASSERT_EQ(positionsVector[0],position2);
    ASSERT_EQ(positionsVector[1],position1);
    ASSERT_EQ(positionsVector[2],position0);
}

TEST_F(APositionsVectorTest, Slice) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    ASSERT_EQ(positionsVector.entity(0).positionsRef(), positions.segment(0,3));
    ASSERT_EQ(positionsVector.entity(1).positionsRef(), positions.segment(3,3));
    ASSERT_EQ(positionsVector.entity(2).positionsRef(), positions.segment(6,3));

    ASSERT_EQ(positionsVector.slice({0,1}).positionsRef(),positionsVector.entity(0).positionsRef());
    ASSERT_EQ(positionsVector.slice({1,1}).positionsRef(),positionsVector.entity(1).positionsRef());
    ASSERT_EQ(positionsVector.slice({2,1}).positionsRef(),positionsVector.entity(2).positionsRef());

    ASSERT_EQ(positionsVector.slice({0,2}).positionsRef(), positions.segment(0,6));
    ASSERT_EQ(positionsVector.slice({1,2}).positionsRef(), positions.segment(3,6));
    ASSERT_EQ(positionsVector.slice({0,3}).positionsRef(), positions.segment(0,9));

    ASSERT_EQ(positionsVector.slice({0,3}).positionsRef(),positionsVector.positionsRef());

    EXPECT_DEATH(positionsVector.entity(-1).positionsRef(),"");
    EXPECT_DEATH(positionsVector.slice({0,4}).positionsRef(),"");
}

TEST_F(APositionsVectorTest, Translate) {
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector positionsVector(positions);

    VectorXd expected[2] = {VectorXd(9),VectorXd(9)};
    expected[0] << 1,2,3,5,6,7,7,8,9;
    expected[1] << 1,2,3,6,7,8,8,9,10;

    positionsVector.entity(1).translate({1,1,1});

    ASSERT_EQ(positionsVector.positionsAsEigenVector(),expected[0]);
    positionsVector.slice({1,2}).translate({1,1,1});

    ASSERT_EQ(positionsVector.positionsAsEigenVector(),expected[1]);
}


TEST_F(APositionsVectorTest, RotateAllCounterClockwise) {
    VectorXd positions(12);
    positions <<\
    0,0,0,\
    1,0,0,\
    0,1,0,\
    0,0,1;
    PositionsVector p(positions);

    double angle = 120.*ConversionFactors::deg2rad; // 120° counterclockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    p.rotateAroundOrigin(angle,axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    0,0,0,\
    0,0,1,\
    1,0,0,\
    0,1,0;

    ASSERT_TRUE(p.positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTest, RotateAllClockwise) {
    VectorXd positions(12);
    positions <<\
    0,0,0,\
    1,0,0,\
    0,1,0,\
    0,0,1;
    PositionsVector p(positions);

    double angle = -120.*ConversionFactors::deg2rad; // 120° clockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    p.rotateAroundOrigin(angle,axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    0,0,0,\
    0,1,0,\
    0,0,1,\
    1,0,0;

    ASSERT_TRUE(p.positionsAsEigenVector().isApprox(expectedPositions));
}


TEST_F(APositionsVectorTest, RotateSliceCounterClockwise) {
    VectorXd positions(12);
    positions <<\
    0,0,0,\
    1,0,0,\
    0,1,0,\
    0,0,1;
    PositionsVector p(positions);

    double angle = 120.*ConversionFactors::deg2rad; // 120° counterclockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    p.slice({2,2}).rotateAroundOrigin(angle,axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    0,0,0,\
    1,0,0,\
    1,0,0,\
    0,1,0;

    ASSERT_TRUE(p.positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTest, RotateSliceClockwise) {
    VectorXd positions(12);
    positions <<\
    0,0,0,\
    1,0,0,\
    0,1,0,\
    0,0,1;
    PositionsVector p(positions);

    double angle = -120.*ConversionFactors::deg2rad; // 120° clockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    p.slice({2,2}).rotateAroundOrigin(angle,axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    0,0,0,\
    1,0,0,\
    0,0,1,\
    1,0,0;

    ASSERT_TRUE(p.positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTest, UglyReturnAndResetHack){
    VectorXd positions(9);
    positions<< position0, position1, position2;
    PositionsVector p(positions);

    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p.position(0),position0);
    ASSERT_EQ(p.slice({1, 2}).position(0),position1);
    ASSERT_EQ(p.position(0),position0);

    ASSERT_EQ(p.position(0),position0);
    ASSERT_EQ(p.slice({1, 2}).position(0, false),position1);
    ASSERT_EQ(p.position(0),position1);
    ASSERT_EQ(p.position(0),position0);

}