//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMPENALIZEDLEASTSQUARESFIT_H
#define RTQC_BSPLINEFROMPENALIZEDLEASTSQUARESFIT_H

#include "BSplineGenerator.h"
#include "BSpline.h"

/* Generates a B-Spline curve with a specified number of control points approximating the data points with loose ends.
 * "Loose ends" means that the first and last control point not need to coincide with the first and last data point of
 * the set.
 * Thus, all n+1 control points are optimzed.
 * Still, the resulting B-Spline curve passes through the first and last control point since the first and last basis
 * functions are equal to unity at the ends of the domain (N(u=0)=1, N(u=1)=1, therefore C(u=0)=P0 and C(u=1)=Plast).
 * Penalization can be turned on by specifying a lambda value > 0 and difference orders can be specified by kappa.
 */

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
