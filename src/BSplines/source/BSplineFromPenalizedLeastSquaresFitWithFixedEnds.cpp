//
// Created by Michael Heuer on 10.05.16.
//

#include <iostream>

#include "BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h"
#include "BSplineTools.h"

// FIXED ENDS: constrained to go through the last data point

BSplineFromPenalizedLeastSquaresFitWithFixedEnds::BSplineFromPenalizedLeastSquaresFitWithFixedEnds(const Eigen::MatrixXd & dataPoints,
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


void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::initializeGenerator(){

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
  generateQVectors();
  calculateConstantTermsMatrix();
  calculateCoefficientMatrix();

  calculateFiniteDifferenceMatrix();

  initializeSolver();
  generateControlPointMatrix();
}

/*! changes the number of polynomial segments of the splines which are equal to the number of control points -1 */
void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::setNumberOfControlPoints(const unsigned numberOfControlPoints) {
  n_=numberOfControlPoints-1;
  initializeGenerator();
}

void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::generateQVectors() {
  Q_.resize(m_+1,dim_);
  Q_.setZero();

  for (unsigned g = 1; g <= m_-1 ; ++g) {
    Q_.row(g) = R_.row(g);
    double rec1= bsBasis_.evaluate(0, p_, n_, U_, uBar_(g));
    double rec2= bsBasis_.evaluate(n_, p_, n_, U_, uBar_(g));
    for (unsigned l = 0; l < dim_; ++l) {
      Q_(g,l)-= rec1*R_(0,l);
      Q_(g,l)-= rec2*R_(m_,l);
    }
  }
}

void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::calculateConstantTermsMatrix() {
  Qmat_.resize(n_-1,dim_);
  for (unsigned i = 1; i <= n_-1; ++i) {
    Eigen::VectorXd accum(dim_,1);
    accum.setZero();

    for (unsigned g = 1; g <= m_-1 ; ++g) {
      accum += bsBasis_.evaluate(i, p_, n_, U_, uBar_(g)) * Q_.row(g);
    }
    Qmat_.row(i-1) = accum;
  }
}

void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::calculateCoefficientMatrix() {
  Nmat_.resize(m_-1,n_-1);
  for (unsigned g = 1; g <= m_-1; ++g) {
    for (unsigned i = 1; i <= n_-1; ++i) {
      Nmat_(g-1,i-1) = bsBasis_.evaluate(i, p_, n_, U_, uBar_(g));
    }
  }
}

void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::initializeSolver() {
  //svdOfNmat_.compute(Nmat_.transpose() * Nmat_, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svdOfNmat_.compute(Nmat_.transpose() * Nmat_ + lambda_ * (DeltaMat_.transpose() * DeltaMat_),
                     Eigen::ComputeThinU | Eigen::ComputeThinV);
}

void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::generateControlPointMatrix() {
  P_.resize(n_+1,dim_);
  P_.setZero();
  Eigen::MatrixXd Pinternal(n_-1,dim_);
  Pinternal=svdOfNmat_.solve(Qmat_);

  // construct control point vector P_
  P_.row(0)  = R_.row(0);
  for (unsigned i = 1; i <= n_-1 ; ++i) { P_.row(i) = Pinternal.row(i-1); };
  P_.row(n_) = R_.row(m_);
}


void BSplineFromPenalizedLeastSquaresFitWithFixedEnds::calculateFiniteDifferenceMatrix() {
  assert(kappa_ <= n_);

  DeltaMat_.resize(n_-1-kappa_,n_-1);//DeltaMat_.resize(n_+1-kappa,n_+1);
  DeltaMat_.setZero();

  for (unsigned i = 0; i < n_-kappa_-2 ; ++i) {
    for (unsigned j = 0; i < n_-2 ; ++j) {
      DeltaMat_(i,j) = BSplineTools::differenceOperator(i,j,kappa_);
    }
  }
}
