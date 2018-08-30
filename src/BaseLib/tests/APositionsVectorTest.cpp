//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include <sstream>
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
        positionX = Vector3d(-1,-1,-1);

        VectorXd tmp(9);
        tmp << position0, position1, position2;
        positions = tmp;

        VectorXd tmp2(6);
        tmp2 << position0, position1;
        smallPositions = tmp2;
    }
    Vector3d position0,position1,position2,positionX;
    VectorXd positions;
    VectorXd smallPositions;
};

TEST_F(APositionsVectorTest, NumberOfEntities){
    PositionsVector p(positions);

    ASSERT_EQ(p.numberOfEntities(),3);
}

TEST_F(APositionsVectorTest, IndexOperator){
    PositionsVector p(positions);

    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p[1],position1);
    ASSERT_EQ(p[2],position2);
    ASSERT_EQ(p[-1],position2);
    ASSERT_EQ(p[-2],position1);
    ASSERT_EQ(p[-3],position0);

    EXPECT_DEATH(p[3],"");
    EXPECT_DEATH(p[-4],"");
}

TEST_F(APositionsVectorTest, Insert){
    PositionsVector p(smallPositions);

    p.insert(positionX,1);

    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p[1],positionX);
    ASSERT_EQ(p[2],position1);
}


TEST_F(APositionsVectorTest, Prepend) {
    PositionsVector p(smallPositions);

    p.prepend(positionX);

    ASSERT_EQ(p[0],positionX);
    ASSERT_EQ(p[1],position0);
    ASSERT_EQ(p[2],position1);
}

TEST_F(APositionsVectorTest, Append) {
    PositionsVector p(smallPositions);

    p.append(positionX);

    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p[1],position1);
    ASSERT_EQ(p[2],positionX);
}

TEST_F(APositionsVectorTest, Print) {
    VectorXd positions(9);
    Vector3d a(0,-2.01238,7.13);
    Vector3d b(1e-10,-9.8987,-12.238);
    Vector3d c(-99.9,99.9,9.99);
    positions<< a, b, c;
    PositionsVector positionsVector(positions);
    std::ostringstream os;
    os << positionsVector;
    ASSERT_EQ(os.str()," 1    0.00000  -2.01238   7.13000\n 2    0.00000  -9.89870 -12.23800\n 3  -99.90000  99.90000   9.99000\n");
}

TEST_F(APositionsVectorTest, PermuteAll) {
    PositionsVector p(positions);

    VectorXi pvec(3);
    pvec << 0,2,1;
    PermutationMatrix<Dynamic> perm(pvec);
    p.permute(perm);

    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p[1],position2);
    ASSERT_EQ(p[2],position1);
}

TEST_F(APositionsVectorTest, PermuteAllCyclic) {
    PositionsVector p(positions);

    VectorXi pvec(3);
    pvec << 2,0,1;
    PermutationMatrix<Dynamic> perm(pvec);
    p.permute(perm);

    ASSERT_EQ(p[0],position1);
    ASSERT_EQ(p[1],position2);
    ASSERT_EQ(p[2],position0);
}

TEST_F(APositionsVectorTest, PermuteSlice) {
    PositionsVector p(positions);

    VectorXi p1(2);
    p1 << 1,0;
    p.slice({1,2}).permute(PermutationMatrix<Dynamic>(p1));
    ASSERT_EQ(p[0],position0);
    ASSERT_EQ(p[1],position2);
    ASSERT_EQ(p[2],position1);

    VectorXi p2(3);
    p2 << 2,0,1;
    p.permute(PermutationMatrix<Dynamic>(p2));
    ASSERT_EQ(p[0],position2);
    ASSERT_EQ(p[1],position1);
    ASSERT_EQ(p[2],position0);
}

TEST_F(APositionsVectorTest, Slice) {
    PositionsVector p(positions);

    ASSERT_EQ(p.entity(0).positionsRef(), positions.segment(0,3));
    ASSERT_EQ(p.entity(1).positionsRef(), positions.segment(3,3));
    ASSERT_EQ(p.entity(2).positionsRef(), positions.segment(6,3));

    ASSERT_EQ(p.slice({0,1}).positionsRef(),p.entity(0).positionsRef());
    ASSERT_EQ(p.slice({1,1}).positionsRef(),p.entity(1).positionsRef());
    ASSERT_EQ(p.slice({2,1}).positionsRef(),p.entity(2).positionsRef());

    ASSERT_EQ(p.slice({0,2}).positionsRef(), positions.segment(0,6));
    ASSERT_EQ(p.slice({1,2}).positionsRef(), positions.segment(3,6));
    ASSERT_EQ(p.slice({0,3}).positionsRef(), positions.segment(0,9));

    ASSERT_EQ(p.slice({0,3}).positionsRef(),p.positionsRef());

    EXPECT_DEATH(p.entity(-1).positionsRef(),"");
    EXPECT_DEATH(p.slice({0,4}).positionsRef(),"");
}

TEST_F(APositionsVectorTest, Translate) {
    PositionsVector p(positions);

    VectorXd expected[2] = {VectorXd(9),VectorXd(9)};
    expected[0] << 1,2,3,5,6,7,7,8,9;
    expected[1] << 1,2,3,6,7,8,8,9,10;

    p.entity(1).translate({1,1,1});
    ASSERT_EQ(p.positionsAsEigenVector(),expected[0]);

    p.slice({1,2}).translate({1,1,1});
    ASSERT_EQ(p.positionsAsEigenVector(),expected[1]);
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
    Vector3d axis = {1,1,1};

    p.rotateAroundOrigin(angle,axis);

    VectorXd expectedPositions(12);
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
    Vector3d axis = {1,1,1};

    p.rotateAroundOrigin(angle,axis);

    VectorXd expectedPositions(12);
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
    Vector3d axis = {1,1,1};

    p.slice({2,2}).rotateAroundOrigin(angle,axis);

    VectorXd expectedPositions(12);
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
    Vector3d axis = {1,1,1};

    p.slice({2,2}).rotateAroundOrigin(angle,axis);

    VectorXd expectedPositions(12);
    expectedPositions << \
    0,0,0,\
    1,0,0,\
    0,0,1,\
    1,0,0;

    ASSERT_TRUE(p.positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTest, ResetAutomatic){
    PositionsVector p(positions);

    ASSERT_EQ(p.position(0),position0);
    ASSERT_EQ(p.slice({1,2}, Reset::Automatic).position(0),position1);
    ASSERT_EQ(p.position(0),position0);
}

TEST_F(APositionsVectorTest, ResetOnFinished){
    PositionsVector p(positions);

    ASSERT_EQ(p.position(0),position0);
    ASSERT_EQ(p.slice({1,2}, Reset::Manual).position(0),position1);
    ASSERT_EQ(p.position(0),position1);
}

TEST_F(APositionsVectorTest, ResetManual){
    PositionsVector p(positions);

    ASSERT_EQ(p.position(0),position0);
    ASSERT_EQ(p.slice({1,2}, Reset::Manual).position(0,Usage::Finished),position1);
    ASSERT_EQ(p.position(0),position1);
    p.resetRef(); // manual reset
    ASSERT_EQ(p.position(0),position0);
}


TEST_F(APositionsVectorTest, OutputStreamOperator){
    PositionsVector p(positions);

    std::stringstream ss;
    ss << p.position(0).transpose()
    << std::endl
    << p.slice({1,2}, Reset::Automatic).positionsRef().transpose()
    << std::endl
    << p.position(0).transpose();

    std::string expected = "1 2 3\n4 5 6 7 8 9\n1 2 3";
    ASSERT_STREQ(expected.c_str(), ss.str().c_str());
}
