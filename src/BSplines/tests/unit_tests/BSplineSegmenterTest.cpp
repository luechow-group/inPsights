//
// Created by Michael Heuer on 01.07.16.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "BSplineSegmenter.h"
#include "BSpline.h"
#include "BSplineFromControlPolygon.h"
#include "ContainerConverter.h"

using namespace testing;

class ABSplineRootSegmenterTest : public Test {
public:
  void SetUp() override { }
};

TEST_F(ABSplineRootSegmenterTest, CorrectSegments) {
  Eigen::Matrix<double,5,1> mat;
  mat << 2, 1, 0, 1, 2;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  BSpline bs= myBSplineGenerator_.generateBSpline();

  std::vector<double> cuts;
  cuts.push_back(0.5);

  BSplineSegmenter segmenter;
  std::vector<BSpline> segments = segmenter.segment(bs, cuts);

  //std::cout << "number of segments "<< segments.size() << std::endl;

  Eigen::VectorXd knotVector0 = segments[0].getKnotVector(0);
  Eigen::VectorXd knotVector1 = segments[1].getKnotVector(0);

  Eigen::VectorXd controlPoints0 = segments[0].getControlPointMatrix(0);
  Eigen::VectorXd controlPoints1 = segments[1].getControlPointMatrix(0);

  //ContainerConverter::printVectorXdForMathematica(knotVector0); std::cout << std::endl;
  //ContainerConverter::printVectorXdForMathematica(controlPoints0); std::cout << std::endl;
  //std::cout << std::endl;
  //ContainerConverter::printVectorXdForMathematica(knotVector1); std::cout << std::endl;
  //ContainerConverter::printVectorXdForMathematica(controlPoints1); std::cout << std::endl;
  //std::cout << std::endl;


  Eigen::VectorXd refKnotVector0;
  refKnotVector0.resize(7);
  refKnotVector0 << 0, 0, 0, 2/3.0, 1, 1, 1;
  ASSERT_TRUE(knotVector0.isApprox(refKnotVector0));

  Eigen::VectorXd refControlPoints0;
  refControlPoints0.resize(4);
  refControlPoints0 << 2, 1, 0.25, 0.25;
  ASSERT_TRUE(controlPoints0.isApprox(refControlPoints0));

  Eigen::VectorXd refKnotVector1;
  refKnotVector1.resize(7);
  refKnotVector1 << 0, 0, 0, 1/3.0, 1, 1, 1;
  ASSERT_TRUE(knotVector1.isApprox(refKnotVector1));

  Eigen::VectorXd refControlPoints1;
  refControlPoints1.resize(4);
  refControlPoints1 << 0.25, 0.25, 1, 2;
  ASSERT_TRUE(controlPoints1.isApprox(refControlPoints1));
}