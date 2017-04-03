//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMLEASTSQUARESFIT_H
#define RTQC_BSPLINEFROMLEASTSQUARESFIT_H

#include "BSplineGenerator.h"
#include "BSpline.h"

/* Generates a B-Spline curve with a specified number of control points approximating the data points with fixed ends.
 * "Fixed ends" means that the first and last control point are constrained to coincide with the first and last data
 * point of
 * the set. (see the NURBS book by Piegl 1997)
 * Thus, all but the first and last control points (n-1 in total) are optimized.
 * Still, the resulting B-Splines curve passes through the first and last control point since the first and last basis
 * functions are equal to unity at the ends of the domain (N(u=0)=1, N(u=1)=1, therefore C(u=0)=P_0=R_0 and C(u=1)
 * =P_last=R_last).
 * Penalization can be turned on by specifying a lambda value > 0 and difference orders can be specified by kappa.
 * */

class BSplineFromPenalizedLeastSquaresFitWithFixedEnds : public BSplineGenerator {
public:
  BSplineFromPenalizedLeastSquaresFitWithFixedEnds(const Eigen::MatrixXd & dataPoints,
                             const unsigned numberOfControlPoints,
                             const unsigned splineDegree = 3,
                             const bool uniformKnotVector = true,
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
