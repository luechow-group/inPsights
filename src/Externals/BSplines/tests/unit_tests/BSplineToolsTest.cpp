//
// Created by Michael Heuer on 17.06.16.
//

#include <gtest/gtest.h>
#include <Eigen/Core>

#include "BSplineTools.h"
#include "BSpline.h"
#include "BSplineFromControlPolygon.h"


using namespace testing;

class ABSplineToolsTest : public Test {
public:
  unsigned p;
  BSpline bs;

  unsigned l;
  unsigned idx;

  void SetUp() override {
  }
};

TEST_F(ABSplineToolsTest, LeftEqualCorrectIdDegree1) {

  p = 1;

  Eigen::MatrixXd data(3,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);
}

TEST_F(ABSplineToolsTest, LeftCorrectIdDegree1) {

  p = 1;

  Eigen::MatrixXd data(3,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfLeftDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfLeftDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);
}

TEST_F(ABSplineToolsTest, RightCorrectIdDegree1) {

  p = 1;

  Eigen::MatrixXd data(3,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfRightDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfRightDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfRightDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfRightDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfRightDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);
}

TEST_F(ABSplineToolsTest, RightOrEqualCorrectIdDegree1) {

  p = 1;

  Eigen::MatrixXd data(3,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,1);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,2);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);
}

TEST_F(ABSplineToolsTest, LeftEqualCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);
  data.row(3) =Eigen::Vector2d(4,0);
  data.row(4) =Eigen::Vector2d(5,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);
}

TEST_F(ABSplineToolsTest, LeftCorrectIdDegree3) {

  p = 3;

  Eigen::MatrixXd data(5,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);
  data.row(3) =Eigen::Vector2d(4,0);
  data.row(4) =Eigen::Vector2d(5,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfLeftDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfLeftDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);
}

TEST_F(ABSplineToolsTest, RightCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);
  data.row(3) =Eigen::Vector2d(4,0);
  data.row(4) =Eigen::Vector2d(5,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfRightDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfRightDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfRightDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);

  l = BSplineTools::findIdxOfRightDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);

  l = BSplineTools::findIdxOfRightDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);
}

TEST_F(ABSplineToolsTest, RightOrEqualCorrectIdDegree3) {

  p = 3;

  Eigen::MatrixXd data(5,2);
  data.row(0) =Eigen::Vector2d(1,1);
  data.row(1) =Eigen::Vector2d(2,0);
  data.row(2) =Eigen::Vector2d(3,1);
  data.row(3) =Eigen::Vector2d(4,0);
  data.row(4) =Eigen::Vector2d(5,1);

  BSplineFromControlPolygon myBSplineGenerator_(data,p,true);
  bs= myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,3);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.4999, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5, p, bs.getKnotVector(0));
  ASSERT_EQ(l,4);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5001, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(1.0, p, bs.getKnotVector(0));
  ASSERT_EQ(l,5);
}
