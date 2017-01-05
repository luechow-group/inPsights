//
// Created by Michael Heuer on 18.06.16.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "BSplineKnotInserter.h"
#include "BSplineTools.h"
#include "BSpline.h"
#include "BSplineFromControlPolygon.h"

#include "ContainerConverter.h"

using namespace testing;

class ABSplineKnotInserterTest : public Test {
public:
  unsigned p;
  BSpline bs;
  BSplineKnotInserter bsKnotInserter;


  void SetUp() override {
  }
};

TEST_F(ABSplineKnotInserterTest, InsertNonExistingKnotDegree3 ) {

  p = 3;

  Eigen::MatrixXd data(5,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);
  data.row(3) =Eigen::Vector2d(4,0);
  data.row(4) =Eigen::Vector2d(5,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();
  bsKnotInserter.insertKnotByReference(0.25, bs);

  Eigen::VectorXd vecResult = bs.getKnotVector(0);
  Eigen::MatrixXd matResult = bs.getControlPointMatrix(0);


  //ContainerConverter::printVectorXdForMathematica(vecResult);
  //ContainerConverter::printMatrixXdForMathematica(matResult);

  Eigen::VectorXd vecCheck(10);
  vecCheck << 0,0,0,0,0.25,0.5,1,1,1,1;
  ASSERT_TRUE(vecResult.isApprox(vecCheck));

  Eigen::MatrixXd matCheck(6,2);
  matCheck << 1,1, 1.5,0.5, 2.25,0.25, 3.25,0.75, 4,0, 5,1;
  ASSERT_TRUE(matResult.isApprox(matCheck));
}
