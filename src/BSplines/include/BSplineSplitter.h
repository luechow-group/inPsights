//
// Created by Michael Heuer on 19.06.16.
//

#ifndef RTQC_BSPLINESPLITTER_H
#define RTQC_BSPLINESPLITTER_H

#include <Eigen/Core>
#include "BSpline.h"
#include "BSplineKnotInserter.h"

/*! Splits a B-spline curve at a give parameter.
 * */
class BSplineSplitter {
public:
  BSplineSplitter();
  std::pair<BSpline, BSpline>
  split(const double u, const BSpline &bs, const std::pair<bool, bool> normalizeKnotVectors);

private:
  BSplineKnotInserter bsKnotInserter_;
};

#endif //RTQC_BSPLINESPLITTER_H
