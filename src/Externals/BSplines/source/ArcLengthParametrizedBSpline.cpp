//
// Created by Michael Heuer on 03.06.16.
//

#include "ArcLengthParametrizedBSpline.h"
#include "BSplineFromPointInterpolation.h"
#include "BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h"

#include "NumericalIntegration.h"

ArcLengthParametrizedBSpline::ArcLengthParametrizedBSpline(const BSpline & interpolationSpline)
  :
  interpolationSpline_(interpolationSpline),
  parametrizationSpline_()
{
  interpolationSpline_.calculateDerivatives(1);
}

ArcLengthParametrizedBSpline::ArcLengthParametrizedBSpline(const Eigen::VectorXd & knotVector,
                                                           const Eigen::MatrixXd & controlPoints,
                                                           const unsigned degree,
                                                           const int highestDerivativeOrder)
  :
  interpolationSpline_(knotVector,controlPoints, degree, highestDerivativeOrder),
  parametrizationSpline_()
{
  interpolationSpline_.calculateDerivatives(1);
}


double ArcLengthParametrizedBSpline::infinitisimalLength(double u) {
  return interpolationSpline_.evaluate(u,1).norm();
}

double ArcLengthParametrizedBSpline::integrateArcLength(const double ua, const double ub){
  assert( (0.0 <= ua) && (ua <= ub) && (ub <= 1.0) );

  //if (fabs(ub-ua) < 1E-15) return 0.0;
  //else return NumericalIntegration::integrate(infinitisimalLength,ua,ub,10000);
  return 0.0;
}

void ArcLengthParametrizedBSpline::integratePolynomialSegments() {

  auto n = unsigned(interpolationSpline_.getControlPointMatrix().rows()-1);

  std::vector<double> polynomialSegmentLengths;

  auto U0 = interpolationSpline_.getKnotVector(0);
  for (unsigned i = 0; i < n; ++i) {
    polynomialSegmentLengths.push_back(integrateArcLength(U0(i),U0(i+1)));
  }

}

void ArcLengthParametrizedBSpline::createInverseFunction(){

}

void ArcLengthParametrizedBSpline::parametrizeByArcLength(const unsigned res, const unsigned controlPointNumber){

  // solve integral for defined resolution over the intervall 0,1

  // TODO ADD INTEGRATION AND ZERO POINT SEARCH WITH BRENT DEKKER

  double ua,ub;
  Eigen::VectorXd arcLengthsIntervals(res);
  Eigen::VectorXd cumulativeArcLengths(res);


  arcLengthsIntervals(0) = 0.0;
  cumulativeArcLengths (0) = 0.0;

  for (unsigned i = 1; i < res; ++i) {
    ua = i-1/(res-1);
    ub = i/(res-1);

    arcLengthsIntervals(i)=integrateArcLength(ua,ub);
    cumulativeArcLengths (i) =  cumulativeArcLengths(i-1) + integrateArcLength(ua,ub);
  }

  totalArcLength_ = cumulativeArcLengths.tail(1)(1);//arcLengthsIntervals.sum();
  Eigen::VectorXd relativeLengths = cumulativeArcLengths.array()/totalArcLength_;

  //TODO MAKE INVERSE FUNC
  BSplineFromPointInterpolation bSplineGenerator(relativeLengths,3);
  parametrizationSpline_ = bSplineGenerator.generateBSpline(0);

}

Eigen::VectorXd ArcLengthParametrizedBSpline::evaluate(const double u, const unsigned k) const {
  assert(isArcLengthParametrized_  && "The spline is not parametrized yet. Run parametrizeByArcLength() "
                                        "or try using deriveAndEvaluate.");
  assert( k >= 0 && "Derivative order has to be greater than zero.");
  assert( (0.0 <= u) && (u <= 1.0) );

  double uArcLength = parametrizationSpline_.evaluate(u,0)(0);
  return BSpline::evaluate(uArcLength, k);
}

Eigen::VectorXd ArcLengthParametrizedBSpline::operator()(const double u, const unsigned int k) const {
  return ArcLengthParametrizedBSpline::evaluate(u,k);
}

Eigen::VectorXd ArcLengthParametrizedBSpline::deriveAndEvaluate(const double u, const unsigned k) {
  if (!isArcLengthParametrized_) parametrizeByArcLength();
  assert( k >= 0 && "Derivative order has to be greater than zero.");
  assert( (0.0 <= u) && (u <= 1.0) );

  double uArcLength = parametrizationSpline_.deriveAndEvaluate(u,0)(0);//TODO - geht das so mit den Ableitungen oder
  // braucht man einen std::vector<BSpline> der alle ableitungen bogenlängenparametrisiert?
  return BSpline::deriveAndEvaluate(uArcLength, k);
}

Eigen::VectorXd ArcLengthParametrizedBSpline::deriveAndEvaluateNaive(const double u, const unsigned int k) {
  if (!isArcLengthParametrized_) parametrizeByArcLength();
  assert( k >= 0 && "Derivative order has to be greater than zero.");
  assert( (0.0 <= u) && (u <= 1.0) );

  double uArcLength = parametrizationSpline_.deriveAndEvaluate(u,0)(0);//TODO - geht das so mit den Ableitungen oder
  // braucht man einen std::vector<BSpline> der alle ableitungen bogenlängenparametrisiert?
  return BSpline::deriveAndEvaluateNaive(uArcLength, k);
}
