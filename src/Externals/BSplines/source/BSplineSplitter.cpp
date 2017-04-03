//
// Created by Michael Heuer on 19.06.16.
//

#include "BSplineSplitter.h"
#include "BSplineTools.h"
#include "ContainerConverter.h"

BSplineSplitter::BSplineSplitter()
  :
   bsKnotInserter_()
{ }


std::pair<BSpline, BSpline>
BSplineSplitter::split(const double u, const BSpline &bs, const std::pair<bool, bool> normalizeKnotVectors) {

  assert( (u > 0.0) && (u < 1.0) );

  unsigned p = bs.getDegree();
  BSpline bsInserted;
  bsInserted = bs;

  unsigned l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(u, p, bs.getKnotVector(0));
  unsigned it = l;
  unsigned multiplicity = 0;

  while(bs.getKnotVector(0)[it] == u) {
    multiplicity++;
    it++;
  }

  assert(multiplicity <= p && "Spline cannot be disconnected");

  for (unsigned i = 0; i < p-multiplicity; ++i) {
    bsKnotInserter_.insertKnotByReference(u, bsInserted);
  }

  if (multiplicity > 0){

    Eigen::VectorXd U(bsInserted.getKnotVector(0));
    Eigen::MatrixXd P(bsInserted.getControlPointMatrix(0));

    // left b-spline curve
    Eigen::VectorXd Uleft(l+p+1);
    Uleft.head(l+p) = U.head(l+p);
    Uleft(l+p) = u;

    if(normalizeKnotVectors.first) BSplineTools::normalizeKnotVector(Uleft);

    Eigen::MatrixXd Pleft(P.topRows(Uleft.size()-(p+2)+1));

    // right b-spline curve
    Eigen::VectorXd Uright(U.size()-l+1); // length of U - the elements that occur before the split idx
    Uright(0) = u;
    Uright.tail(U.size()-(l)) = U.tail(U.size()-(l));

    if(normalizeKnotVectors.second) BSplineTools::normalizeKnotVector(Uright);

    Eigen::MatrixXd Pright(P.bottomRows(Uright.size()-(p+2)+1));

    return std::make_pair(BSpline(Uleft,Pleft,p,0),BSpline(Uright,Pright,p,0));
  }
  else {

    Eigen::VectorXd U(bsInserted.getKnotVector(0));
    Eigen::MatrixXd P(bsInserted.getControlPointMatrix(0));

    // left b-spline curve
    Eigen::VectorXd Uleft(l+p+1+1);
    Uleft.head(l+p+1) = U.head(l+p+1);
    Uleft(l+p+1) = u;

    if(normalizeKnotVectors.first) BSplineTools::normalizeKnotVector(Uleft);

    //Eigen::MatrixXd Pleft(P.topRows(l+1));
    Eigen::MatrixXd Pleft(P.topRows(Uleft.size()-(p+2)+1));

    // right b-spline curve
    Eigen::VectorXd Uright(U.size()-(l+1)+1); // length of U - the elements that occur before the split idx
    Uright(0) = u;
    Uright.tail(U.size()-(l+1)) = U.tail(U.size()-(l+1));

    if(normalizeKnotVectors.second) BSplineTools::normalizeKnotVector(Uright);

    //Eigen::MatrixXd Pright(P.bottomRows(P.rows()-l));
    Eigen::MatrixXd Pright(P.bottomRows(Uright.size()-(p+2)+1));

    return std::make_pair(BSpline(Uleft,Pleft,p,0),BSpline(Uright,Pright,p,0));
  }

}

