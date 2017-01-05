//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMLEASTSQUARESFIT_H
#define RTQC_BSPLINEFROMLEASTSQUARESFIT_H

#include "BSplineGenerator.h"
#include "BSpline.h"

/*! Generates a B-spline curve with a specified number of control points which approximates the data points.
 *  (see the NURBS book by Piegl 1997)
 * */
// fixed ends
class BSplineFromPenalizedLeastSquaresFitWithFixedEnds : public BSplineGenerator {
public:
  BSplineFromPenalizedLeastSquaresFitWithFixedEnds(const Eigen::MatrixXd & dataPoints,
                             const unsigned numberOfControlPoints,
                             const unsigned splineDegree = 3,
                             const bool uniformKnotVector = false,
                             const double lambda = 0,
                             const unsigned kappa = 2);

  void setNumberOfControlPoints(const unsigned numberOfControlPoints);

private:
  void initializeGenerator() override;

  /*! eq. (9.63) in the NURBS book */
  void generateQVectors();

  /*! eq. (9.67) in the NURBS book */
  void calculateConstantTermsMatrix();

  /*! eq. (9.66) in the NURBS book */
  void calculateCoefficientMatrix();

  /*! Calculates the control point by least squares curve approximation
   * according to eq. (9.65) in the NURBS book
   * The values of P0=Q0 and Pn=Qn are pre- and appended to the control point vector */
  void generateControlPointMatrix() override;

  void calculateFiniteDifferenceMatrix();

  void initializeSolver();

  bool uniformKnotVector_;
  double lambda_;
  unsigned kappa_;
  Eigen::MatrixXd Q_,Qmat_,Nmat_,DeltaMat_;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
};


#endif //RTQC_BSPLINEFROMLEASTSQUARESFIT_H
