//
// Created by Michael Heuer on 30.06.16.
//

#ifndef RTQC_BSPLINEEXTREMAFINDER_H
#define RTQC_BSPLINEEXTREMAFINDER_H

#include "BSplineRootFinder.h"
#include "BSpline.h"

/*! Analyses B-spline curves and allows to detect motifs (minima, maxima, saddle-points)
 * - The treatment of roots and root intervals in the first derivative at the B-spline curve ends is arbitrary and
 * needs to  be revised, especially, with respect to an application in the path-processing scheme of a continuously
 * growing exploration path (buffering of path segments etc.).
 * - A pure reorientation means that the entire B-spline curve is completely flat so that the first derivative is
 * always zero. Accordingly, the B-spline curve contains a root interval ranging over u=[0,1].
 * //TODO In this case a bool isCompletelyFlat would be maybe more clear.
 * */
class BSplineStationaryPointFinder{
public:
  BSplineStationaryPointFinder(const BSpline& bs);

  std::vector<double> getMinima(const unsigned dimIdx, const double firstDerivativeZeroThreshold = 1e-13);
  std::vector<double> getMaxima(const unsigned dimIdx, const double firstDerivativeZeroThreshold = 1e-13);
  std::vector<double> getIncSaddlePoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold = 1e-13);
  std::vector<double> getDecSaddlePoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold = 1e-13);

  std::vector<std::pair<double,double> > getBasins(const unsigned dimIdx, const double firstDerivativeZeroThreshold =
  1e-13);
  std::vector<std::pair<double,double> > getPlateaus(const unsigned dimIdx, const double firstDerivativeZeroThreshold
  = 1e-13);
  std::vector<std::pair<double,double> > getIncSaddles(const unsigned dimIdx, const double
  firstDerivativeZeroThreshold = 1e-13);
  std::vector<std::pair<double,double> > getDecSaddles(const unsigned dimIdx, const double
  firstDerivativeZeroThreshold = 1e-13);
  std::vector<std::pair<double,double> > getPureReorientation(const unsigned dimIdx, const double
  firstDerivativeZeroThreshold = 1e-13);

private:

  void findStationaryPoints(const unsigned dimIdx, const double firstDerivativeZeroThreshold = 1e-13);
  void clear();

  unsigned p_,currentDim_;
  double firstDerivativeZeroThreshold_;
  BSpline bsOrigCopy_;

  std::vector<double> singleRootsOfFirstDervivative_,repeatedRootsOfFirstDerivative_, minima_,maxima_,incSaddlePoints_,
    decSaddlePoints_;
  std::vector<std::pair<double,double > > rootIntervalsOfFirstDerivative_,basins_,plateaus_,incSaddles_, decSaddles_,
    pureReorientations_;
};

#endif //RTQC_BSPLINEEXTREMAFINDER_H
