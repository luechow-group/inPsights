//
// Created by Michael Heuer on 01.07.16.
//

#include "BSplineSegmenter.h"
#include "BSpline.h"

BSplineSegmenter::BSplineSegmenter()
  : bsSplitter_() {}

std::vector<BSpline> BSplineSegmenter::segment(const BSpline &bs, const std::vector<double> &splitParams){

  assert(std::is_sorted(splitParams.begin(), splitParams.end()) && "The split parameters have to be in ascending order.");
  assert(*std::min_element(splitParams.begin(), splitParams.end()) >= 1e-14);
  assert(*std::max_element(splitParams.begin(), splitParams.end()) <= 1-1e-14);

  std::vector<BSpline> segments;
  std::pair<BSpline,BSpline> splitResult;
  splitResult.second = bs;

  if(splitParams.size() >=1 ) {

    // iterate and stop before the last cut parameter
    for (unsigned i = 0; i < splitParams.size() - 1; ++i) {
      //don't normalize the knot vector of right spline yet
      splitResult = bsSplitter_.split(splitParams[i], splitResult.second, {true, false});
      segments.push_back(splitResult.first);
    }
    // do last cut and also normalize the knot vector of the right spline
    splitResult = bsSplitter_.split(splitParams.back(), splitResult.second, {true, true});

    segments.push_back(splitResult.first);
    segments.push_back(splitResult.second);

  }
  return segments;
}