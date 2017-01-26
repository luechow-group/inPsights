//
// Created by Michael Heuer on 26.01.17.
//


#include <gmock/gmock.h>
#include <Eigen/Core>

#include "Optimization/LBFGSOptimizer.h"

using namespace testing;
using namespace Eigen;

class AHessianFromGradientCalculatorTest : public Test {
public:
  void SetUp() override {

  }
};

TEST_F(AHessianFromGradientCalculatorTest, OriginCentered2DParabola){

  Eigen::Vector2d samplePosition(1,1);
  Eigen::Vector2d exactGradient(2, 2);
  Eigen::Matrix2d exacthessian(2,0, 0,2);

  // calculate hessian from gradient


  ASSERT_TRUE(false);
}