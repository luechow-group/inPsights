//
// Created by Michael Heuer on 17.06.16.
//

#ifndef RTQC_BSPLINEROOTFINDER_H
#define RTQC_BSPLINEROOTFINDER_H

#include <Eigen/Core>

#include "BSpline.h"
#include "BSplineKnotInserter.h"
#include "BSplineSplitter.h"


struct Interval{
  explicit Interval(){
    empty = true;
  }
  explicit Interval(unsigned startIdx)
  {
    start = startIdx;
    empty = false;
  }
  explicit Interval(unsigned startIdx, unsigned endIdx) {
    start = startIdx;
    end= endIdx;
    empty = false;
  }

  unsigned length() const {
    return end - start + 1;
  }

  unsigned start;
  unsigned end;
  bool empty;
};

/*! Finds roots and root intervals in selected dimension of a B-spline curve (that must have a polynomial degree >= 2).
 * */
class BSplineRootFinder {
public:
  BSplineRootFinder(const BSpline& bs);

  void findRoots(const unsigned dimIdx, const double zeroThreshold = 1e-13);

  std::vector<double> getSingleRoots(){ return singleRoots_;}

  std::vector<double> getRepeatedRoots(){ return repeatedRoots_;}

  std::vector<std::pair<double,double> > getRootIntervals(){ return rootIntervals_;}

  std::vector<BSpline> getSegmentedBSplines(){ return segmentedBSplines_; }

private:

  void findIndicesOfControlPointsEqualToZero(const unsigned dimIdx);
  void findIntervalsOfControlPointsEqualToZero();

  void findObviousRoots();
  void findObviousRootsOrRootIntervalsAtStart();
  void findObviousRootsOrRootIntervalsAtEnd();
  void findObviousRootsOrRootIntervalsInMiddle(const Interval &zeroInterval);

  void findAllUnobviousRootsInBSplineSegment(BSpline &bsSeg, const unsigned dimIdx);
  std::pair<bool,double> findLeftmostUnobviousRootInBSplineSegment(BSpline &bs, const unsigned dimIdx);
  double calculateRootApproximant(const unsigned i, const unsigned dimIdx) const;
  bool checkConvergence(const unsigned i);

  unsigned p_,numberOfControlPointsInOriginalBSpline_,maxIterations_;
  double convergenceThreshold_,zeroThreshold_;
  Eigen::VectorXd uRootEstimates_;


  const BSpline bsOrig_;
  BSpline bsTemp_;

  Eigen::VectorXd Uref_;
  Eigen::MatrixXd Pref_;

  BSplineKnotInserter bsKnotInserter_;
  BSplineSplitter bsSplitter_;

  std::vector<unsigned> indicesOfControlPointsEqualToZero_;
  std::vector<Interval> limitingIndicesOfIntervalsOfControlPointsEqualToZero_;

  std::vector<double> singleRoots_,repeatedRoots_;
  std::vector<std::pair<double,double> > rootIntervals_;

  std::vector<BSpline> segmentedBSplines_;

  bool firstControlPointIsZero_,lastControlPointIsZero_;
  bool splineStartsWithZeroControlPoint_;

};

#endif //RTQC_BSPLINEROOTFINDER_H
