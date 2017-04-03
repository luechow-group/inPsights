//
// Created by Michael Heuer on 19.05.16.
//

#ifndef RTQC_BSPLINEMERGER_H
#define RTQC_BSPLINEMERGER_H

#include <Eigen/Core>
#include <Eigen/SVD>
#include <vector>

#include "BSpline.h"
#include "ContainerConverter.h"

/* Merging of two p-th order B-Splines via lagrange optimization and knot adjustment:
 *  |Tai, C.-L., Hu, S.-M. & Huang, Q.-X.
 *  |Approximate merging of B-spline curves via knot adjustment and constrained optimization.
 *  |Comput. Des. 35, 893–899 (2003).
 * doi: 10.1016/S0010-4485(02)00176-8
 * However, the paper is written very unclear. I suggest reading my masters thesis.
 */

class BSplineMerger {
public:
  BSplineMerger();

  /* merges two splines */
  BSpline merge(const BSpline &bs1, const BSpline &bs2);

  /* merges two splines with keeping the last control point of the first spline fixed*/
  BSpline connectToBSpline1(const BSpline &bs1, const BSpline &bs2);

  /* merges two splines with keeping the first control point of the second spline fixed*/
  BSpline connectToBSpline2(const BSpline &bs1, const BSpline &bs2);

  void printResultsForMathematica();

private:
  BSpline mergeWithConstraints(const BSpline &bs1, const std::vector<double> &constrainedParams1,
                               const BSpline &bs2, const std::vector<double> &constrainedParams2);

  bool kroneckerDelta(const unsigned i, const unsigned j) const;
  double prefactor(const unsigned idx, const unsigned zerothOrderIdx, const unsigned k,
                   const Eigen::MatrixXd &zerothOrderControlPointMatrix, const Eigen::VectorXd &zerothOrderKnotVector) const;

  void initialize(const BSpline &bs1,const BSpline &bs2);
  Eigen::MatrixXd calculateKV();
  Eigen::MatrixXd calculateKW();
  Eigen::MatrixXd calculateGV();
  Eigen::MatrixXd calculateHW();
  Eigen::MatrixXd calculateIV();
  Eigen::MatrixXd calculateJW();
  Eigen::MatrixXd calculateIpV();
  Eigen::MatrixXd calculateJppW();
  Eigen::MatrixXd constructNmat();
  Eigen::MatrixXd calculateKconst();
  Eigen::MatrixXd constructConstMat();
  void initializeSolver();
  void solveLinearEquationSystem();
  void generateShiftedControlPointVectorsOfSpline1and2();
  void generateControlPointVectorOfMergedSpline();

  void adjustParametrization();

  Eigen::VectorXd reverseKnotVector(Eigen::VectorXd U);

  Eigen::VectorXd adjustKnotVector(const Eigen::VectorXd &Uleft,
                                   const unsigned nLeft,
                                   const Eigen::VectorXd &Uright);
  void adjustKnotVectorsOfBothSplines();

  void createMergedKnotVector();
  Eigen::VectorXd generateDerivativeControlPoint(const unsigned idx,
                                                 const unsigned k,
                                                 const Eigen::MatrixXd &zerothOrderControlPointMatrix,
                                                 const Eigen::VectorXd &zerothOrderKnotVector) const;



  Eigen::MatrixXd adjustShiftedControlPoints(const Eigen::MatrixXd &shiftedControlPoints,
                                             const Eigen::VectorXd &originalKnotVector,
                                             const Eigen::VectorXd &adjustedKnotVector,
                                             const unsigned n);
  void adjustShiftedControlPointsOfBothSplines();

  unsigned p_,dim_,m_,o_,mg_,oh_,
    numberOfConstraintsA_,
    numberOfConstraintsB_;
  std::vector<double> constrainedParamsA_, constrainedParamsB_;
  Eigen::VectorXd Vrev_,Wrev_,Vadj_,Wadj_,Wadjrev_,Umerged_;
  Eigen::MatrixXd S_,T_,shifts_,Sshifted_,Tshifted_,Sadj_,Tadj_,Tadjrev_,Pmerged_;

  std::vector<Eigen::VectorXd> Vk_,Wk_;
  std::vector<Eigen::MatrixXd> Sk_,Tk_;

  //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> svdOfNmat_;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
  BSplineBasis bsBasis_;
};

#endif //RTQC_BSPLINEMERGER_H
