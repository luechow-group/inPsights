//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEFROMCONTROLPOLYGON_H
#define RTQC_BSPLINEFROMCONTROLPOLYGON_H

#include "BSplineGenerator.h"

/*! Generates a B-spline curve from a set of data points that are used as control points.
 * */
class BSplineFromControlPolygon : public BSplineGenerator {
public:
  BSplineFromControlPolygon(const Eigen::MatrixXd & dataPoints,
                            const unsigned degree = 3,
                            const bool uniformKnotVector = false);

private:
  bool uniformKnotVector_;

  void initializeGenerator() override;
  void generateControlPointMatrix() override;
};

#endif //RTQC_BSPLINEFROMCONTROLPOLYGON_H
