//
// Created by Michael Heuer on 30.06.16.
//

#include "BSplineStationaryPointFinder.h"
#include "BSplineTools.h"

BSplineStationaryPointFinder::BSplineStationaryPointFinder(const BSpline& bs)
  :
  p_(bs.getDegree()),
  currentDim_(std::numeric_limits<unsigned >::max()),
  bsOrigCopy_(bs)
{
  assert(bs.getDegree() >= 2 && "BSpline must be at least of degree p=2.");
  //bsOrig_.calculateDerivatives(1);
}

void BSplineStationaryPointFinder::clear() {
  singleRootsOfFirstDervivative_.clear();
  minima_.clear();
  maxima_.clear();

  repeatedRootsOfFirstDerivative_.clear();
  decSaddlePoints_.clear();
  incSaddlePoints_.clear();

  rootIntervalsOfFirstDerivative_.clear();
  basins_.clear();
  plateaus_.clear();
  incSaddles_.clear();
  decSaddles_.clear();
  pureReorientations_.clear();
}

void BSplineStationaryPointFinder::findStationaryPoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold){
  clear();

  assert(dimIdx < bsOrigCopy_.getDim());

  currentDim_ = dimIdx;
  firstDerivativeZeroThreshold_=firstDerivativeZeroThreshold;
  BSplineRootFinder bsRootFinderForFirstDerivative(bsOrigCopy_.getDerivativeBSpline(1));

  bsRootFinderForFirstDerivative.findRoots(dimIdx,firstDerivativeZeroThreshold_);
  singleRootsOfFirstDervivative_  = bsRootFinderForFirstDerivative.getSingleRoots();
  repeatedRootsOfFirstDerivative_ = bsRootFinderForFirstDerivative.getRepeatedRoots();
  rootIntervalsOfFirstDerivative_ = bsRootFinderForFirstDerivative.getRootIntervals();

  const auto& Uref = bsOrigCopy_.getKnotVector(0);

  for (auto root : singleRootsOfFirstDervivative_) {

    // determine minima and maxima by sign shift
    unsigned leftKnotIdx = BSplineTools::findIdxOfLeftDomainKnot(root, p_, Uref);
    unsigned rightKnotIdx = BSplineTools::findIdxOfRightDomainKnot(root, p_, Uref);
    double valueOfFirstDerivativeAtLeftKnot = bsOrigCopy_.evaluate(Uref(leftKnotIdx),1)(dimIdx);
    double valueOfFirstDerivativeAtRightKnot = bsOrigCopy_.evaluate(Uref(rightKnotIdx),1)(dimIdx);

    if (root == 0.0) {
      if      (valueOfFirstDerivativeAtRightKnot > 0) minima_.push_back(root);
      else if (valueOfFirstDerivativeAtRightKnot < 0) maxima_.push_back(root);
    }
    if (root == 1.0) {
      if      (valueOfFirstDerivativeAtLeftKnot < 0) minima_.push_back(root);
      else if (valueOfFirstDerivativeAtLeftKnot > 0) maxima_.push_back(root);
    }
    else {
      double valueAtRoot = bsOrigCopy_.deriveAndEvaluate(root, 2)(dimIdx);
      if      (valueAtRoot > 0) minima_.push_back(root);
      else if (valueAtRoot < 0) maxima_.push_back(root);
        // if this point is reached, the value at the root is == 0
      else if ( valueOfFirstDerivativeAtLeftKnot * valueOfFirstDerivativeAtRightKnot > 0) {
        if ((valueOfFirstDerivativeAtLeftKnot > 0) && (valueOfFirstDerivativeAtRightKnot > 0))
          incSaddlePoints_.push_back(root);
        else if ((valueOfFirstDerivativeAtLeftKnot < 0) && (valueOfFirstDerivativeAtRightKnot < 0))
          decSaddlePoints_.push_back(root);
      }
    }
  }

  for (auto root : repeatedRootsOfFirstDerivative_) {

    /* repeated roots at the head and tail are immediately saved as single roots since saddle points cannot exist at
     * the head and tail.
     */

    // use sign change criterion as a sufficient condition to determine if minimum, maximum or saddle-point
    unsigned leftKnotIdx = BSplineTools::findIdxOfLeftDomainKnot(root, p_, Uref);
    unsigned rightKnotIdx = BSplineTools::findIdxOfRightDomainKnot(root, p_, Uref);
    double valueOfFirstDerivativeAtLeftKnot = bsOrigCopy_.evaluate(Uref(leftKnotIdx),1)(dimIdx);
    double valueOfFirstDerivativeAtRightKnot = bsOrigCopy_.evaluate(Uref(rightKnotIdx),1)(dimIdx);

    if ( (valueOfFirstDerivativeAtLeftKnot < 0) && (valueOfFirstDerivativeAtRightKnot > 0) )
      minima_.push_back(root);
    else if ( (valueOfFirstDerivativeAtLeftKnot > 0) && (valueOfFirstDerivativeAtRightKnot < 0) )
      maxima_.push_back(root);
    else if ( valueOfFirstDerivativeAtLeftKnot * valueOfFirstDerivativeAtRightKnot > 0) {
      if ((valueOfFirstDerivativeAtLeftKnot > 0) && (valueOfFirstDerivativeAtRightKnot > 0))
        incSaddlePoints_.push_back(root);
      else if ((valueOfFirstDerivativeAtLeftKnot < 0) && (valueOfFirstDerivativeAtRightKnot < 0))
        decSaddlePoints_.push_back(root);
    }
    else
      throw ("Error: One of the neighboring knots is zero and therefore the current root is part of a root interval. "
        "This should not have been happened.");
  }

  std::sort(minima_.begin(),minima_.end());
  std::sort(maxima_.begin(),maxima_.end());


  for (auto& rootInterval : rootIntervalsOfFirstDerivative_) {

    unsigned leftKnotIdx = BSplineTools::findIdxOfLeftDomainKnot(rootInterval.first, p_, Uref);
    unsigned rightKnotIdx = BSplineTools::findIdxOfRightDomainKnot(rootInterval.second, p_, Uref);
    double valueOfFirstDerivativeAtLeftKnot = bsOrigCopy_.evaluate(Uref(leftKnotIdx),1)(dimIdx);
    double valueOfFirstDerivativeAtRightKnot = bsOrigCopy_.evaluate(Uref(rightKnotIdx),1)(dimIdx);

    // root interval spans the whole spline
    if ((rootInterval.first == 0.0) && (rootInterval.second == 1.0)) {
      pureReorientations_.push_back(rootInterval);
    }
      // root interval at the head
    else if (rootInterval.first == 0.0) {
      if (valueOfFirstDerivativeAtRightKnot > 0) basins_.push_back(rootInterval);
      else if (valueOfFirstDerivativeAtRightKnot < 0) plateaus_.push_back(rootInterval);
      pureReorientations_.push_back(rootInterval);
    }
      // root interval at the tail
    else if (rootInterval.second == 1.0) {
      if (valueOfFirstDerivativeAtLeftKnot < 0) basins_.push_back(rootInterval);
      else if (valueOfFirstDerivativeAtLeftKnot > 0) plateaus_.push_back(rootInterval);
      pureReorientations_.push_back(rootInterval);
    }
      // interior root intervals
    else {
      if ( (valueOfFirstDerivativeAtLeftKnot < 0) && (valueOfFirstDerivativeAtRightKnot > 0) )
        basins_.push_back(rootInterval);
      else if ( (valueOfFirstDerivativeAtLeftKnot > 0) && (valueOfFirstDerivativeAtRightKnot < 0) )
        plateaus_.push_back(rootInterval);
      else if ( valueOfFirstDerivativeAtLeftKnot * valueOfFirstDerivativeAtRightKnot > 0) {
        if ((valueOfFirstDerivativeAtLeftKnot > 0) && (valueOfFirstDerivativeAtRightKnot > 0))
          incSaddles_.push_back(rootInterval);
        else if ((valueOfFirstDerivativeAtLeftKnot < 0) && (valueOfFirstDerivativeAtRightKnot < 0))
          decSaddles_.push_back(rootInterval);
      }
      else pureReorientations_.push_back(rootInterval);
    }
  }
}

std::vector<double> BSplineStationaryPointFinder::getMinima(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
  if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return minima_;
}

std::vector<double> BSplineStationaryPointFinder::getMaxima(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return maxima_;
}

std::vector<double> BSplineStationaryPointFinder::getIncSaddlePoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return incSaddlePoints_;
}

std::vector<double> BSplineStationaryPointFinder::getDecSaddlePoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return decSaddlePoints_;
}

std::vector<std::pair<double,double>> BSplineStationaryPointFinder::getBasins(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return basins_;
}

std::vector<std::pair<double,double>> BSplineStationaryPointFinder::getPlateaus(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return plateaus_;
}

std::vector<std::pair<double,double>> BSplineStationaryPointFinder::getIncSaddles(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return incSaddles_;
}

std::vector<std::pair<double,double>> BSplineStationaryPointFinder::getDecSaddles(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return decSaddles_;
}


std::vector<std::pair<double,double>> BSplineStationaryPointFinder::getPureReorientation(const unsigned dimIdx, const double firstDerivativeZeroThreshold) {
    if ((firstDerivativeZeroThreshold != firstDerivativeZeroThreshold_) || (dimIdx != currentDim_)) findStationaryPoints(dimIdx, firstDerivativeZeroThreshold);
  return pureReorientations_;
}

