//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINE_H
#define RTQC_BSPLINE_H

#include <Eigen/Core>
#include <vector>
#include "BSplineBasis.h"

/*! Contains the standard implementation of a basis-spline curve.
 * Evaluation and getter methods are defined so that knot vectors and control points of derivatives are calculated
 * when they are needed.
 * Defines virtual functions for the ArcLengthParametrizedBSpline class that is unfinished yet.
 * */
class BSpline {
public:
  BSpline();
  virtual ~BSpline(){}

  BSpline(const Eigen::VectorXd& knotVector,
          const Eigen::MatrixXd& controlPoints,
          const unsigned degree,
          int highestDerivativeToCalculate = -1);

  virtual Eigen::VectorXd evaluate(const double u, const unsigned derivativeOrder = 0) const;

  virtual Eigen::VectorXd operator()(const double u, const unsigned derivativeOrder=0) const;

  virtual Eigen::VectorXd deriveAndEvaluate(const double u, const unsigned derivativeOrder = 0);

  virtual Eigen::VectorXd deriveAndEvaluateNaive(const double u, const unsigned derivativeOrder = 0);

  
  void reverseSpline();

  unsigned getDegree() const {
    return p_;
  }

  unsigned getOrder() const {
    return p_+1;
  }

  unsigned getDim() const {
    return dim_;
  }

  unsigned getControlPointNumber() const {
    return unsigned(getControlPointMatrix(0).rows());
  }

  BSpline getDerivativeBSpline(const unsigned derivative);
  
  const Eigen::VectorXd& getKnotVector(const unsigned derivative=0) const;
  
  const std::vector<Eigen::VectorXd>& getKnotVectors() const;
  
  const Eigen::MatrixXd& getControlPointMatrix(const unsigned derivativeOrder = 0) const;
  
  const std::vector<Eigen::MatrixXd>& getControlPointMatrices() const;

  void calculateDerivatives(const unsigned highestDerivativeToCalculate);

  const Eigen::VectorXd& deriveAndGetKnotVector(const unsigned derivative);
  
  const Eigen::MatrixXd& deriveAndGetControlPointMatrix(const unsigned derivativeOrder);

  bool isClampedAndNormed() const;


private:
  Eigen::VectorXd deBoorAlgorithm(const double u,const unsigned i,const unsigned p,const unsigned k=0) const;

  void deriveKnotVectors(const unsigned highestDerivativeToCalculate);
  
  Eigen::VectorXd deriveControlPoint(const unsigned i, const unsigned highestDerivativeToCalculate) const;

  void deriveControlPointMatrices(const unsigned int highestDerivativeToCalculate);

  unsigned findIdxOfLowerOrEqualDomainKnot(const double u, const unsigned k) const;

  unsigned p_,dim_,n_,highestCalculatedDerivative_;
  std::vector<Eigen::VectorXd> Uk_;
  std::vector<Eigen::MatrixXd> Pk_;
  BSplineBasis bsBasis_;
  bool isReversed_;
};

#endif //RTQC_BSPLINE_H
