//
// Created by Michael Heuer on 17.06.16.
//

#include <cmath>

#include "BSplineRootFinder.h"
#include "BSplineTools.h"
#include "ContainerConverter.h"




BSplineRootFinder::BSplineRootFinder(const BSpline& bs)
  :
  p_(bs.getDegree()),
  numberOfControlPointsInOriginalBSpline_((unsigned) bs.getControlPointMatrix(0).rows()-1),
  maxIterations_(10000),
  convergenceThreshold_(1e-13),
  bsOrig_(bs),
  bsTemp_(bs),
  Uref_(bs.getKnotVector(0)),
  Pref_(bs.getControlPointMatrix(0)),
  bsKnotInserter_(),
  bsSplitter_()
{
  assert(bs.isClampedAndNormed() );
}

void BSplineRootFinder::findRoots(const unsigned dimIdx, const double zeroThreshold){
  assert(dimIdx < bsOrig_.getDim());
  zeroThreshold_ = std::abs(zeroThreshold);

  findIndicesOfControlPointsEqualToZero(dimIdx);
  findIntervalsOfControlPointsEqualToZero();
  findObviousRoots();

  for (auto bsSeg : segmentedBSplines_) findAllUnobviousRootsInBSplineSegment(bsSeg, dimIdx);

  std::sort(singleRoots_.begin(), singleRoots_.end());
}

void BSplineRootFinder::findIndicesOfControlPointsEqualToZero(const unsigned dimIdx){
  for (unsigned i = 0; i < Pref_.rows(); ++i) {
    if(std::abs(Pref_(i,dimIdx)) < zeroThreshold_ ) indicesOfControlPointsEqualToZero_.push_back(i);
  }
}

void BSplineRootFinder::findIntervalsOfControlPointsEqualToZero(){
  if(indicesOfControlPointsEqualToZero_.size() > 0) {
    for (unsigned startIdx = 0; startIdx <= indicesOfControlPointsEqualToZero_.size() - 1; ++startIdx) {
      unsigned endIdx = startIdx;
      while ((indicesOfControlPointsEqualToZero_[endIdx] + 1 == indicesOfControlPointsEqualToZero_[endIdx + 1])
             && (endIdx <= indicesOfControlPointsEqualToZero_.size()))
        ++endIdx;

      limitingIndicesOfIntervalsOfControlPointsEqualToZero_.push_back(Interval(indicesOfControlPointsEqualToZero_[startIdx], indicesOfControlPointsEqualToZero_[endIdx]));

      startIdx = endIdx;
    }
  }
}

void BSplineRootFinder::findObviousRoots(){
  /* if no multiple roots are present in the BSpline */
  if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.empty()) segmentedBSplines_.push_back(bsOrig_);

    /* if all control points are equal to zero */
  else if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.front().end == numberOfControlPointsInOriginalBSpline_) rootIntervals_.push_back(std::make_pair(0.0, 1.0));

    /* special root cases are present */
  else {

    // check start and end

    unsigned treatStart = 0;
    unsigned treatEnd = (unsigned) limitingIndicesOfIntervalsOfControlPointsEqualToZero_.size()-1;

    if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.front().start == 0){
      firstControlPointIsZero_ = true;
      treatStart++;
    } else  firstControlPointIsZero_ = false;

    if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.back().end == numberOfControlPointsInOriginalBSpline_) {
      lastControlPointIsZero_ = true;
      treatEnd--;
    } else  lastControlPointIsZero_ = false;


    // analyze the polygon
    if (firstControlPointIsZero_) findObviousRootsOrRootIntervalsAtStart();
    for (unsigned i = treatStart; i <= treatEnd ; ++i) {
      findObviousRootsOrRootIntervalsInMiddle(limitingIndicesOfIntervalsOfControlPointsEqualToZero_[i]);
    }
    if (lastControlPointIsZero_) findObviousRootsOrRootIntervalsAtEnd();
    else segmentedBSplines_.push_back(bsTemp_);
  }
}

void BSplineRootFinder::findObviousRootsOrRootIntervalsAtStart(){
  /* first zero interval is a single root */
  /* double root cannot exist at the start because the third derivative cannot be non-zero */
  if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.front().length() <= p_) singleRoots_.push_back(0.0);
    /* first zero interval is a root interval */
  else {
    unsigned startKnotIdx = p_ + limitingIndicesOfIntervalsOfControlPointsEqualToZero_.front().start;
    unsigned endKnotIdx = p_ + limitingIndicesOfIntervalsOfControlPointsEqualToZero_.front().end - (p_ - 1);

    double uRootIntervalStart = Uref_(startKnotIdx);
    double uRootIntervalEnd = Uref_(endKnotIdx);

    rootIntervals_.push_back(std::make_pair(uRootIntervalStart, uRootIntervalEnd));

    /* split at interval end without rescaling the knot vectors */
    auto bsLeftRightPair = bsSplitter_.split(uRootIntervalEnd, bsTemp_, {false, false});
    bsTemp_ = bsLeftRightPair.second;
  }
}

void BSplineRootFinder::findObviousRootsOrRootIntervalsAtEnd(){
  /* last zero interval is a single root (double root cannot exist at the end per definition) */
  if (limitingIndicesOfIntervalsOfControlPointsEqualToZero_.back().length() <= p_) {
    singleRoots_.push_back(1.0);
    segmentedBSplines_.push_back(bsTemp_);
  }
    /* last zero interval is a root interval */
  else {
    unsigned startKnotIdx = p_ + limitingIndicesOfIntervalsOfControlPointsEqualToZero_.back().start;
    unsigned endKnotIdx = p_ + limitingIndicesOfIntervalsOfControlPointsEqualToZero_.back().end - (p_ - 1);

    double uRootIntervalStart = Uref_(startKnotIdx);
    double uRootIntervalEnd = Uref_(endKnotIdx);

    rootIntervals_.push_back(std::make_pair(uRootIntervalStart, uRootIntervalEnd));

    /* split at interval start without rescaling the knot vectors */
    auto bsLeftRightPair = bsSplitter_.split(uRootIntervalStart, bsTemp_, {false, false});
    segmentedBSplines_.push_back(bsLeftRightPair.first);
  }
}

void BSplineRootFinder::findObviousRootsOrRootIntervalsInMiddle(const Interval &zeroInterval){
  unsigned startKnotIdx = p_ + zeroInterval.start;
  unsigned endKnotIdx = p_ + zeroInterval.end - (p_ - 1);

  /* single roots are treated by the findAllUnobviousRootsInBSplineSegment algorithm */

  /* treat interior repeated root */
  if (zeroInterval.length() == p_) {

    double uRepeatedRoot = Uref_(startKnotIdx);
    repeatedRoots_.push_back(uRepeatedRoot);

    /* split at double root without rescaling the knot vectors */
    auto bsLeftRightPair = bsSplitter_.split(uRepeatedRoot, bsTemp_, {false, false});
    segmentedBSplines_.push_back(bsLeftRightPair.first);
    bsTemp_ = bsLeftRightPair.second;
  }
    /* treat interior root interval */
  else if (zeroInterval.length() > p_) {
    double uRootIntervalStart = Uref_(startKnotIdx);
    double uRootIntervalEnd = Uref_(endKnotIdx);

    rootIntervals_.push_back(std::make_pair(uRootIntervalStart, uRootIntervalEnd));

    // split two times at the start and end of the root interval without rescaling the knot vectors (and without saving
    // the middle spline segements but the outer ones)
    auto bsLeftRightPair1 = bsSplitter_.split(uRootIntervalStart, bsTemp_, {false, false});
    segmentedBSplines_.push_back(bsLeftRightPair1.first);

    auto bsLeftRightPair2 = bsSplitter_.split(uRootIntervalEnd, bsLeftRightPair1.second, {false, false});
    bsTemp_ = bsLeftRightPair2.second;
  }
}

void BSplineRootFinder::findAllUnobviousRootsInBSplineSegment(BSpline &bs, const unsigned dimIdx){

  splineStartsWithZeroControlPoint_ = false;
  double root;
  bool rootFound;
  std::pair<BSpline,BSpline> splitResult;
  splitResult.second = bs;

  // check if first control point is zero and do first step
  if (std::abs(bs.getControlPointMatrix(0)(0,dimIdx)) < zeroThreshold_) splineStartsWithZeroControlPoint_ = true;

  auto result = findLeftmostUnobviousRootInBSplineSegment(bs, dimIdx);
  rootFound   = result.first;
  root        = result.second;

  while( rootFound ){
    singleRoots_.push_back(root);
    splitResult = bsSplitter_.split(root, splitResult.second, {false, false});
    result = findLeftmostUnobviousRootInBSplineSegment(splitResult.second, dimIdx);
    rootFound = result.first;
    root = result.second;
  }
}

std::pair<bool,double> BSplineRootFinder::findLeftmostUnobviousRootInBSplineSegment(BSpline &bs, const unsigned dimIdx) {
  Uref_ = bs.getKnotVector();
  Pref_ = bs.getControlPointMatrix();

  BSpline searchSpline_ = bs;

  unsigned i = 1;
  unsigned iLast = (unsigned) Pref_.rows()-1;

  /* skip leading and trailing zeros */
  if (splineStartsWithZeroControlPoint_){
    while ( fabs(Pref_(i,dimIdx))< convergenceThreshold_ && (i < iLast) ) ++i;
    ++i;
  }
  while ( fabs(Pref_(iLast,dimIdx))< convergenceThreshold_ && (iLast >= 1) ) --iLast;


  uRootEstimates_.resize(0);

  for (unsigned it = 0; it < maxIterations_ ; ++it) {
    double uRootEstimate;

    // find left index of the control point interval in which the sign change occurs
    while ( (i <= iLast) && (Pref_(i-1,dimIdx)*Pref_(i,dimIdx)>0 ) ) ++i;

    if (i > iLast) return std::make_pair(false,-1);

    uRootEstimate = calculateRootApproximant(i, dimIdx);
    uRootEstimates_.conservativeResize(uRootEstimates_.size()+1);
    uRootEstimates_.tail(1)(0) = uRootEstimate;


    if (checkConvergence(i)) {
      if (!splineStartsWithZeroControlPoint_) {splineStartsWithZeroControlPoint_ = true;}
      return std::make_pair(true,uRootEstimate);
    }
    else {
      // if not converged insert knot
      bsKnotInserter_.insertKnotByReference(uRootEstimate, searchSpline_);
      iLast++; // iLast grows because of the knot insertion
      Uref_ = searchSpline_.getKnotVector(0);
      Pref_ = searchSpline_.getControlPointMatrix(0);
    }
  }
  return std::make_pair(false,-1);
}

double BSplineRootFinder::calculateRootApproximant(const unsigned i, const unsigned dimIdx) const{
  double u = BSplineTools::knotAverage(i,p_,Uref_);
  u -= Pref_(i,dimIdx)*(Uref_(i+p_)-Uref_(i))/(Pref_(i,dimIdx)-Pref_(i-1,dimIdx))/(p_);
  return u;
}

bool BSplineRootFinder::checkConvergence(const unsigned i) {

  if (uRootEstimates_.size() < p_) return false;

  Eigen::VectorXd lastKnots = uRootEstimates_.tail(p_);
  double uListMax = lastKnots.maxCoeff();
  double uListMin = lastKnots.minCoeff();

  return (uListMax-uListMin) / std::max(Uref_(i),Uref_(i+p_)) < convergenceThreshold_;
}


