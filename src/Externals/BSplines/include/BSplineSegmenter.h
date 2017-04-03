//
// Created by Michael Heuer on 01.07.16.
//

#ifndef RTQC_BSPLINESEGMENTER_H
#define RTQC_BSPLINESEGMENTER_H

#include "BSplineSplitter.h"

class BSpline;

/*! Convenience class for splitting a B-spline curve multiple times into several segments.
 * */
class BSplineSegmenter{
public:
  BSplineSegmenter();
  std::vector<BSpline> segment(const BSpline &bs, const std::vector<double> &splitParams);

private:
  BSplineSplitter bsSplitter_;
};

#endif //RTQC_BSPLINESEGMENTER_H
