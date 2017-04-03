//
// Created by Michael Heuer on 10.05.16.
//

#include <iostream>

#include "BSplineFromPenalizedLeastSquaresFitWithLooseEnds.h"
#include "BSplineTools.h"


BSplineFromPenalizedLeastSquaresFitWithLooseEnds::BSplineFromPenalizedLeastSquaresFitWithLooseEnds(const Eigen::MatrixXd & dataPoints,
                                                                         const unsigned numberOfControlPoints,
                                                                         const unsigned splineDegree,
                                                                         const bool uniformKnotVector,
                                                                         const double lambda,
                                                                         const unsigned kappa)
  :
  BSplineGenerator(dataPoints, numberOfControlPoints, splineDegree),
  uniformKnotVector_(uniformKnotVector),
  lambda_(lambda),
  kappa_(kappa)
{
  initializeGenerator();
}


void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::initializeGenerator(){

  if(!uniformKnotVector_) {
    //generateParametersByChordLengthMethod();
    generateParametersByCentripetalMethod();
    
    //generateKnotVectorByKnotAveraging();
    generateKnotVectorByDeBoorsMethod();
  }
  else {
    generateParametersByEquallySpacedMethod();
    generateKnotVectorByUniformMethod();
  }
  //generateQVectors();
  //calculateConstantTermsMatrix();
  calculateCoefficientMatrix();
  calculateFiniteDifferenceMatrix();
  initializeSolver();
  generateControlPointMatrix();


}

/*! changes the number of polynomial segments of the splines which are equal to the number of control points -1 */
void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::setNumberOfControlPoints(const unsigned numberOfControlPoints) {
  n_=numberOfControlPoints-1;
  initializeGenerator();
}


void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::calculateCoefficientMatrix() {
  Nmat_.resize(m_+1,n_+1);
  for (unsigned g = 0; g <= m_; ++g) {
    for (unsigned i = 0; i <= n_; ++i) {
      Nmat_(g,i) = bsBasis_.evaluate(i, p_, n_, U_, uBar_(g));
    }
  }
}

void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::initializeSolver() {
  svdOfNmat_.compute(Nmat_.transpose() * Nmat_ + lambda_ * (DeltaMat_.transpose() * DeltaMat_),
                     Eigen::ComputeThinU | Eigen::ComputeThinV);
}

void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::generateControlPointMatrix() {
  P_.resize(n_+1,dim_);
  P_.setZero();
  P_=svdOfNmat_.solve(Nmat_.transpose() * R_);

}


void BSplineFromPenalizedLeastSquaresFitWithLooseEnds::calculateFiniteDifferenceMatrix() {
  assert(kappa_ <= n_);

  DeltaMat_.resize(n_+1-kappa_,n_+1);
  DeltaMat_.setZero();

  for (unsigned i = 0; i <= n_-kappa_ ; ++i) {
    for (unsigned j = 0; j <= n_ ; ++j) {
      DeltaMat_(i,j) = BSplineTools::differenceOperator(i,j,kappa_);
    }
  }
}
