//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMPENALIZEDLEASTSQUARESFIT_H
#define RTQC_BSPLINEFROMPENALIZEDLEASTSQUARESFIT_H

#include "BSplineGenerator.h"
#include "BSpline.h"

/*! Generates a B-spline curve with a specified number of control points which approximates the data points.
 *  (see the NURBS book by Piegl 1997)
 * */
// loose ends
class BSplineFromPenalizedLeastSquaresFitWithLooseEnds : public BSplineGenerator {
public:
  BSplineFromPenalizedLeastSquaresFitWithLooseEnds(const Eigen::MatrixXd & dataPoints,
                                      const unsigned numberOfControlPoints,
                                      const unsigned splineDegree = 3,
                                      const bool uniformKnotVector = true,
                                      const double lambda = 0,
                                      const unsigned kappa = 2);

  void setNumberOfControlPoints(const unsigned numberOfControlPoints);

private:
  void initializeGenerator() override;


  /*! eq. (9.66) in the NURBS book */
  void calculateCoefficientMatrix();

  /*! Calculates the control point by least squares curve approximation
   * DIFFERENTLY */
  void generateControlPointMatrix() override;

  void calculateFiniteDifferenceMatrix();


  void initializeSolver();

  bool uniformKnotVector_;
  double lambda_;
  unsigned kappa_;
  Eigen::MatrixXd Nmat_,DeltaMat_;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
};


#endif //RTQC_BSPLINEFROMPENALIZEDLEASTSQUARESFIT_H
