//
// Created by Michael Heuer on 17.06.16.
//

#include "BSplineKnotInserter.h"
#include "BSpline.h"
#include "BSplineTools.h"

#include "iostream"

BSplineKnotInserter::BSplineKnotInserter(){}

/*!
 * Boehms algorithm for knot insertion
 * the knot must be >0 and <1 */
void BSplineKnotInserter::insertKnotByReference(const double uInsert, BSpline &bs) const {
  assert( (uInsert > 0.0) && (uInsert < 1.0) );

  unsigned p = bs.getDegree();
  unsigned dim = bs.getDim();

  Eigen::VectorXd Uold = bs.getKnotVector(0);
  Eigen::MatrixXd Pold = bs.getControlPointMatrix(0);
  unsigned numberOfControlPoints = (unsigned) Pold.rows();
  unsigned numberOfKnots = (unsigned) Uold.size();

  // find knot span
  unsigned l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(uInsert, p, Uold);

  // compute new control points. Only the control points l-p+1 to l change
  Eigen::MatrixXd Pnew(numberOfControlPoints+1,dim);
  Pnew.topRows(l-p+1) = Pold.topRows(l-p+1);
  Pnew.bottomRows(numberOfControlPoints-l) = Pold.bottomRows(numberOfControlPoints-l);

  double alpha;
  for (unsigned i = l-p+1; i <= l; ++i) {
    alpha = (uInsert-Uold(i))/(Uold(i+p)-Uold(i));
    Pnew.row(i) = (1-alpha)*Pold.row(i-1) + alpha*Pold.row(i);
  }

  // insert u into knot vector
  Eigen::VectorXd Unew(numberOfKnots+1);
  Unew.head(l+1) = Uold.head(l+1);
  Unew(l+1) = uInsert;
  Unew.tail(numberOfKnots-l-1) = Uold.tail(numberOfKnots-l-1);

  bs = BSpline(Unew,Pnew,p,0);
}

BSpline BSplineKnotInserter::insertKnotByCopy(const double u, const BSpline &bs) const {
  auto bsCopy = bs;
  BSplineKnotInserter::insertKnotByReference(u, bsCopy);
  return bsCopy;
}
