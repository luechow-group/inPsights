//
// Created by Michael Heuer on 14.06.16.
//

#include "gmock/gmock.h"
#include "BSpline.h"
#include "BSplineMerger.h"
#include "BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h"
#include "BSplineFromControlPolygon.h"

using namespace testing;

class ABSplineMergerTest : public Test {
public:
  unsigned p;
  BSpline bsMerged;
  void SetUp() override {

    p = 3;
    Eigen::MatrixXd data1(6,2);
    data1.row(0) =Eigen::Vector2d(-4.0,-4.00);
    data1.row(1) =Eigen::Vector2d(-3.0,-2.00);
    data1.row(2) =Eigen::Vector2d(-2.0,-0.75);
    data1.row(3) =Eigen::Vector2d(-1.0,-0.25);
    data1.row(4) =Eigen::Vector2d(-0.5,-0.05);
    data1.row(5) =Eigen::Vector2d(+0.0,-0.00);


    Eigen::MatrixXd data2(6,2);
    data2.row(0) =Eigen::Vector2d(+0.0,+0.10);
    data2.row(1) =Eigen::Vector2d(+0.5,+0.05);
    data2.row(2) =Eigen::Vector2d(+1.0,+0.25);
    data2.row(3) =Eigen::Vector2d(+2.0,+0.75);
    data2.row(4) =Eigen::Vector2d(+3.0,+2.00);
    data2.row(5) =Eigen::Vector2d(+4.0,+4.00);

    BSpline bs1= BSplineFromControlPolygon(data1,p,true).generateBSpline();
    BSpline bs2= BSplineFromControlPolygon(data2,p,true).generateBSpline();

    BSplineMerger bsMerger_;
    bsMerged = bsMerger_.connectToBSpline1(bs1, bs2);

  }
};

TEST_F(ABSplineMergerTest,CheckConnectToFirst) {
ASSERT_THAT(bsMerged.evaluate(0.5)(0), DoubleNear(0,10E-15));
ASSERT_THAT(bsMerged.evaluate(0.5)(1), DoubleNear(0,10E-15));
}
