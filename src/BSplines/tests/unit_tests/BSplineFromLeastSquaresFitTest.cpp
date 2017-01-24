//
// Created by Michael Heuer on 17.10.16.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "BSpline.h"
#include "BSplineFromPenalizedLeastSquaresFitWithLooseEnds.h"
#include "BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h"
#include "ContainerConverter.h"

using namespace testing;

class ABSplineFromLeastSquaresFitTest : public Test {
public:
  unsigned p;
  BSpline bs;
  void SetUp() override {
  }
};

TEST_F(ABSplineFromLeastSquaresFitTest, Penalized ) {

  p = 3;

  Eigen::MatrixXd data(10,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,2);
  data.row(2) =Eigen::Vector2d(3,3);
  data.row(3) =Eigen::Vector2d(4,4);
  data.row(4) =Eigen::Vector2d(5,5);
  data.row(5) =Eigen::Vector2d(6,6);
  data.row(6) =Eigen::Vector2d(7,7);
  data.row(7) =Eigen::Vector2d(8,8);
  data.row(8) =Eigen::Vector2d(9,9);
  data.row(9) =Eigen::Vector2d(10,10);

  //BSplineFromPenalizedLeastSquaresFitWithFixedEnds myBSplineGenerator_(data,10,p,true,0,2);
  BSplineFromPenalizedLeastSquaresFitWithLooseEnds myBSplineGenerator_(data,5,p,true,0.5,2);
  bs= myBSplineGenerator_.generateBSpline();

  Eigen::VectorXd vecResult = bs.getKnotVector(0);
  Eigen::MatrixXd matResult = bs.getControlPointMatrix(0);


  ContainerConverter::printVectorXdForMathematica(vecResult);
  ContainerConverter::printMatrixXdForMathematica(matResult);
  //std::cout << std::endl;

  Eigen::VectorXd vecCheck(10);
  vecCheck << 0,0,0,0,0.25,0.5,1,1,1,1;
  //ASSERT_TRUE(vecResult.isApprox(vecCheck));

  Eigen::MatrixXd matCheck(10,2);
  matCheck << 1,1, 1.5,0.5, 2.25,0.25, 3.25,0.75, 4,0, 5,1;
  //ASSERT_TRUE(matResult.isApprox(matCheck));

  assert(false && "Test missing");

}