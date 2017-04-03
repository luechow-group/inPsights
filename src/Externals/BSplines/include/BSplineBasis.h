//
// Created by Michael Heuer on 10.05.16.
//

#ifndef RTQC_BSPLINEBASIS_H
#define RTQC_BSPLINEBASIS_H

#include <Eigen/Core>

/*! Evaluate the B-Spline basis functions that are defined by the knot vector
 * */
class BSplineBasis{

public:
  BSplineBasis(){};

  /* Evaluates the basis function with the Cox-DeBoor-Mansfield Recurrence Relation
   * - i can range from 0 to n with n being the number of polynomial segements spline
   * - the number of control points or number of basis functions = n+1
   * - n is the number of polynomial segments
   * - the knot vector can be of any degree k
   * - u ranges form 0 to 1 */

  double evaluate(const unsigned i,
                  const unsigned pk,
                  const unsigned numberOfSplineSegments,
                  const Eigen::VectorXd &knotvector,
                  const double u) const;
};

#endif //RTQC_BSPLINEBASIS_H
