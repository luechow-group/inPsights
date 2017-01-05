//
// Created by Michael Heuer on 17.06.16.
//

#ifndef RTQC_BSPLINETOOLS_H
#define RTQC_BSPLINETOOLS_H

#include <Eigen/Core>

/*! Contains methods that can be used by several B-spline classes.
 * */
class BSplineTools{
public:
  BSplineTools();

  static unsigned findIdxOfLeftOrEqualDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U);
  static unsigned findIdxOfLeftDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U);
  static unsigned findIdxOfRightDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U);
  static unsigned findIdxOfRightOrEqualDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U);


  static double knotAverage(const unsigned i, const unsigned p, const Eigen::VectorXd &U);

  static void normalizeKnotVector(Eigen::VectorXd & knotVector);

  static Eigen::VectorXd normalizedKnotVector(const Eigen::VectorXd &knotVector);

  static double rescaledKnot(const double knot,
                            const std::pair<double,double> oldLim,
                            const std::pair<double,double> newLim);

  static void rescaleKnotVector(Eigen::VectorXd &knotVector,
                                const std::pair<double,double> oldLim,
                                const std::pair<double,double> newLim);

  static Eigen::VectorXd rescaledKnotVector(const Eigen::VectorXd &knotVector,
                                           const std::pair<double,double> oldLim,
                                           const std::pair<double,double> newLim);

  /* finite difference matrix for penalizing BSplines */
  static int differenceOperator(unsigned i, unsigned j, unsigned kappa);
};

#endif //RTQC_BSPLINETOOLS_H
