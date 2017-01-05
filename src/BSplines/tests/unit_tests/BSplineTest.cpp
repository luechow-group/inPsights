//
// Created by Michael Heuer on 16.06.16.
//

#include <gmock/gmock.h>
#include <Eigen/Core>

#include "BSpline.h"
#include "BSplineFromControlPolygon.h"
#include "ContainerConverter.h"

using namespace testing;

class ABSplineTest : public Test {
public:
  unsigned p;
  unsigned dim;
  BSpline bs;
  void SetUp() override {

    p = 3;
    Eigen::MatrixXd data(5,2);
    data.row(0) =Eigen::Vector2d(1,1);
    data.row(1) =Eigen::Vector2d(1.5,-0.5);
    data.row(2) =Eigen::Vector2d(3,1);
    data.row(3) =Eigen::Vector2d(4.5,0.5);
    data.row(4) =Eigen::Vector2d(5,-1);

    BSplineFromControlPolygon bsGenerator(data,p, true);
    bs = bsGenerator.generateBSpline(-1);

    //ContainerConverter::printVectorXdForMathematica(bs.getKnotVector(1));
    //ContainerConverter::printMatrixXdForMathematica(bs.getControlPointMatrix(1));
  }
};

TEST_F(ABSplineTest,CorrectDerivatives0) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(1., 1.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 0).isApprox(ref));
  ref = Eigen::Vector2d(2.552, 0.312);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 0).isApprox(ref));
  ref = Eigen::Vector2d(3.448, 0.584);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 0).isApprox(ref));
  ref = Eigen::Vector2d(5., -1.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 0).isApprox(ref));
}

TEST_F(ABSplineTest,CorrectDerivatives1) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(3., -9.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 1).isApprox(ref));
  ref = Eigen::Vector2d(4.44, 2.04);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 1).isApprox(ref));
  ref = Eigen::Vector2d(4.44, 0.12);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 1).isApprox(ref));
  ref = Eigen::Vector2d(3., -9.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 1).isApprox(ref));
}

TEST_F(ABSplineTest,CorrectDerivatives2) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(6., 54.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 2).isApprox(ref));
  ref = Eigen::Vector2d(1.2, 1.2);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 2).isApprox(ref));
  ref = Eigen::Vector2d(-1.2, -15.6);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 2).isApprox(ref));
  ref = Eigen::Vector2d(-6., -30.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 2).isApprox(ref));
}

TEST_F(ABSplineTest,CorrectDerivatives3) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(-12., -132.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -132.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -36.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -36.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 3).isApprox(ref));
}

TEST_F(ABSplineTest,CorrectDerivatives4) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 4).isApprox(ref));
}

TEST_F(ABSplineTest,CorrectDerivatives5) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.0, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.4, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(0.6, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.deriveAndEvaluate(1.0, 5).isApprox(ref));
}



TEST_F(ABSplineTest,CompareEvaluationMethods) {
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.deriveAndEvaluateNaive(0.0, k).isApprox(bs.deriveAndEvaluate(0.0, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.deriveAndEvaluateNaive(1/3.0, k).isApprox(bs.deriveAndEvaluate(1/3.0, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.deriveAndEvaluateNaive(0.5, k).isApprox(bs.deriveAndEvaluate(0.5, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.deriveAndEvaluateNaive(0.6, k).isApprox(bs.deriveAndEvaluate(0.6, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.deriveAndEvaluateNaive(1.0, k).isApprox(bs.deriveAndEvaluate(1.0, k)));
  }
}

