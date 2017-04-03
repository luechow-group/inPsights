//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMPOINTINTERPOLATION_H
#define RTQC_BSPLINEFROMPOINTINTERPOLATION_H

#include "BSplineGenerator.h"

/*! Generates a B-spline curve interpolating the data points.
 * (see the NURBS book by Piegl 1997)
 * */
class BSplineFromPointInterpolation : public BSplineGenerator {
public:
  BSplineFromPointInterpolation(const Eigen::MatrixXd & dataPoints,
                                const unsigned degree = 3,
                                const bool uniformKnotVector = false);

private:
  void initializeGenerator() override;
  void generateControlPointMatrix() override;
  void calculateCoefficientMatrix();
  void initializeSolver();

  bool uniformKnotVector_;
  Eigen::MatrixXd Nmat_;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
};

#endif //RTQC_BSPLINEFROMPOINTINTERPOLATION_H
