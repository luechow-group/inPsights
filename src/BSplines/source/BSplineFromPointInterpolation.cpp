//
// Created by Michael Heuer on 10.05.16.
//

#include "BSplineFromPointInterpolation.h"

#include <iostream>

BSplineFromPointInterpolation::BSplineFromPointInterpolation(const Eigen::MatrixXd & dataPoints,
                                                             const unsigned degree,
                                                             const bool uniformKnotVector)
  : BSplineGenerator(dataPoints, degree), uniformKnotVector_(uniformKnotVector)
{
  initializeGenerator();
}

void BSplineFromPointInterpolation::initializeGenerator(){
  if (uniformKnotVector_){
    generateKnotVectorByUniformMethod();
  }
  else {
    //generateParametersByChordLengthMethod();
    generateParametersByCentripetalMethod();
    generateKnotVectorByKnotAveraging();
  }
  calculateCoefficientMatrix();
  initializeSolver();
  generateControlPointMatrix();
}

void BSplineFromPointInterpolation::calculateCoefficientMatrix() {
  Nmat_.resize(n_+1,n_+1);
  for (unsigned g = 0; g <= m_; ++g) {
    for (unsigned i = 0; i <= n_; ++i) {
      Nmat_(g,i) = bsBasis_.evaluate(i, p_, n_, U_, uBar_(g));
    }
  }
}

void BSplineFromPointInterpolation::initializeSolver() {
  //svdOfNmat_.compute(Nmat_);
  svdOfNmat_.compute(Nmat_, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

void BSplineFromPointInterpolation::generateControlPointMatrix() {
  P_.resize(n_+1,dim_);
  P_.setZero();

  P_ = svdOfNmat_.solve(R_);
}