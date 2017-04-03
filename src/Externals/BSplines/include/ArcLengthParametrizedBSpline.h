//
// Created by Michael Heuer on 03.06.16.
//

#ifndef RTQC_ARCLENGTHPARAMETRIZEDBSPLINE_H
#define RTQC_ARCLENGTHPARAMETRIZEDBSPLINE_H

#include "BSpline.h"

class ArcLengthParametrizedBSpline : public BSpline {

public:
  ArcLengthParametrizedBSpline(const BSpline & interpolationSpline);
  ArcLengthParametrizedBSpline(const Eigen::VectorXd & knotVector,
                               const Eigen::MatrixXd & controlPoints,
                               const unsigned degree,
                               const int highestDerivativeOrder = -1);

  double infinitisimalLength(double u);

  void integratePolynomialSegments();

  Eigen::VectorXd evaluate(const double u, const unsigned derivativeOrder = 0) const override;

  Eigen::VectorXd operator()(const double u, const unsigned derivativeOrder=0) const override;

  Eigen::VectorXd deriveAndEvaluate(const double u, const unsigned derivativeOrder = 0) override;

  Eigen::VectorXd deriveAndEvaluateNaive(const double u, const unsigned derivativeOrder = 0) override;


  double totalArcLength_;

private:
  BSpline interpolationSpline_, parametrizationSpline_;

  double integrateArcLength(const double ua, const double ub);
  void createInverseFunction();

  void parametrizeByArcLength(const unsigned resolution = 21, const unsigned controlPointNumber =10);
  bool isArcLengthParametrized_;
};

#endif //RTQC_ARCLENGTHPARAMETRIZEDBSPLINE_H
