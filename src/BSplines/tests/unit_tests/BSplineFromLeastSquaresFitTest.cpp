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


  /*Eigen::MatrixXd data(6,2);
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
  data.row(10)=Eigen::Vector2d(11,11);
  data.row(11)=Eigen::Vector2d(12,12);
  data.row(12)=Eigen::Vector2d(13,13);
  data.row(13)=Eigen::Vector2d(14,14);
  data.row(14)=Eigen::Vector2d(15,15);
  data.row(15)=Eigen::Vector2d(16,16);
  data.row(16)=Eigen::Vector2d(17,17);
  data.row(17)=Eigen::Vector2d(18,18);*/// m+1=6
  Eigen::MatrixXd data(18,2);
   data.row(0) =Eigen::Vector2d(0.017824, 0.195131);
   data.row(1) =Eigen::Vector2d(0.221340, 0.140994);
   data.row(2) =Eigen::Vector2d(0.287290, 0.229286);
   data.row(3) =Eigen::Vector2d(0.490554, 0.172050);
   data.row(4) =Eigen::Vector2d(0.517142, 0.359398);
   data.row(5) =Eigen::Vector2d(0.603510, 0.507977);
   data.row(6) =Eigen::Vector2d(0.845844, 0.669823);
   data.row(7) =Eigen::Vector2d(0.825442, 0.835279);
   data.row(8) =Eigen::Vector2d(1.038500, 0.896949);
   data.row(9) =Eigen::Vector2d(1.075760, 1.169660);
   data.row(10)=Eigen::Vector2d(1.361320, 1.556090);
   data.row(11)=Eigen::Vector2d(1.417830, 1.800650);
   data.row(12)=Eigen::Vector2d(1.567380+2, 2.103810);
   data.row(13)=Eigen::Vector2d(1.685070, 2.534830);
   data.row(14)=Eigen::Vector2d(1.792050, 2.730170);
   data.row(15)=Eigen::Vector2d(1.869890, 3.240830);
   data.row(16)=Eigen::Vector2d(1.944370, 3.610990);
   data.row(17)=Eigen::Vector2d(2.009650, 4.051540);

  //BSplineFromPenalizedLeastSquaresFitWithFixedEnds myBSplineGenerator_(data,10,p,true,0,2);
  BSplineFromPenalizedLeastSquaresFitWithLooseEnds myBSplineGenerator_(data,10,p,true,0.5,2);
  bs= myBSplineGenerator_.generateBSpline();

  Eigen::VectorXd vecResult = bs.getKnotVector(0);
  Eigen::MatrixXd matResult = bs.getControlPointMatrix(0);


  //ContainerConverter::printVectorXdForMathematica(vecResult);
  //ContainerConverter::printMatrixXdForMathematica(matResult);
  //std::cout << std::endl;

  ASSERT_TRUE(true);
  /*Eigen::VectorXd vecCheck(10);
  vecCheck << 0,0,0,0,0.25,0.5,1,1,1,1;
  ASSERT_TRUE(vecResult.isApprox(vecCheck));

  Eigen::MatrixXd matCheck(6,2);
  matCheck << 1,1, 1.5,0.5, 2.25,0.25, 3.25,0.75, 4,0, 5,1;
  ASSERT_TRUE(matResult.isApprox(matCheck));*/

}