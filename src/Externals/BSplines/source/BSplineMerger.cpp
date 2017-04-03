//
// Created by Michael Heuer on 19.05.16.
//

#include "BSplineMerger.h"
#include <iostream>

BSplineMerger::BSplineMerger() { }

void BSplineMerger::initialize(const BSpline & bs1, const BSpline & bs2){
  assert( bs1.getDegree() == bs2.getDegree() && "Both B-Splines must have the same degree.");
  assert( bs1.getDim() == bs2.getDim() && "Both B-Splines must have the same dimension." );
  assert(bs1.isClampedAndNormed() && bs2.isClampedAndNormed() && "Both B-Splines must be normed");

  p_ = bs1.getDegree();
  dim_ = bs1.getDim();

  // spline 1
  Vk_= bs1.getKnotVectors();
  Sk_= bs1.getControlPointMatrices();
  S_ = Sk_[0];
  m_ = (unsigned) S_.rows()-1;

  // spline 2
  Wk_= bs2.getKnotVectors();
  Tk_= bs2.getControlPointMatrices();
  T_ = Tk_[0];
  o_ = (unsigned) T_.rows()-1;
}


BSpline BSplineMerger::merge(const BSpline &bs1, const BSpline &bs2){
  return mergeWithConstraints(bs1, {}, bs2, {});
}

BSpline BSplineMerger::connectToBSpline1(const BSpline &bs1, const BSpline &bs2){
  return mergeWithConstraints(bs1, {1.0}, bs2, {});
}

BSpline BSplineMerger::connectToBSpline2(const BSpline &bs1, const BSpline &bs2){
  return mergeWithConstraints(bs1, {}, bs2, {0});
}

BSpline BSplineMerger::mergeWithConstraints(const BSpline &bs1,
                                            const std::vector<double> &constrainedParams1,
                                            const BSpline &bs2,
                                            const std::vector<double> &constrainedParams2) {
  initialize(bs1,bs2);
  constrainedParamsA_ = constrainedParams1;
  constrainedParamsB_ = constrainedParams2;
  numberOfConstraintsA_ = (unsigned) constrainedParamsA_.size();
  numberOfConstraintsB_ = (unsigned) constrainedParamsB_.size();

  assert( numberOfConstraintsA_ + numberOfConstraintsB_ < p_ &&
           "The total number of constrained points must be less than the spline degree. "
             "Otherwise no solution for the linear system of equations exists.");

  initializeSolver();

  // shift splines so that they fulfill the precise merging conditions
  solveLinearEquationSystem();
  generateShiftedControlPointVectorsOfSpline1and2();

  // adjust and merge knot vectors
  adjustKnotVectorsOfBothSplines();
  createMergedKnotVector();

  // adjust control points
  adjustShiftedControlPointsOfBothSplines();

  // do actual merge
  generateControlPointVectorOfMergedSpline();

  //printResultsForMathematica();

  BSpline bs12 (Umerged_,Pmerged_,p_,-1);
  return bs12;
}


bool BSplineMerger::kroneckerDelta(const unsigned i, const unsigned j) const {
  return i == j;
}

double BSplineMerger::prefactor(const unsigned idx,
                                const unsigned zerothOrderIdx, const unsigned k,
                                const Eigen::MatrixXd &zerothOrderControlPointMatrix,
                                const Eigen::VectorXd &zerothOrderKnotVector) const {
  long n = zerothOrderControlPointMatrix.rows()-1;

  if( (idx<=n-k) ) {
    if (k == 0) {
      return kroneckerDelta(idx,zerothOrderIdx)?1.0:0.0;
    } else if (zerothOrderKnotVector(idx+p_+1) == zerothOrderKnotVector(idx+k)) {
      return 0;
    } else {
      return double(p_+1-k)/(zerothOrderKnotVector(idx+p_+1)-zerothOrderKnotVector(idx+k))
             * (prefactor(idx + 1, zerothOrderIdx, k - 1, zerothOrderControlPointMatrix, zerothOrderKnotVector)
                - prefactor(idx, zerothOrderIdx, k - 1, zerothOrderControlPointMatrix, zerothOrderKnotVector));
    }
  } else {
    return 0;
  }
}

Eigen::VectorXd BSplineMerger::generateDerivativeControlPoint(const unsigned idx,
                                                              const unsigned k,
                                                              const Eigen::MatrixXd & zerothOrderControlPointMatrix,
                                                              const Eigen::VectorXd & zerothOrderKnotVector) const {
  Eigen::VectorXd controlPoint = Eigen::VectorXd::Zero(dim_);

  assert( (idx <= zerothOrderControlPointMatrix.rows()-1-k) && "Index out of bounds: only n-k control points exist");
  for (unsigned zeroOrderIdx =0 ; zeroOrderIdx < zerothOrderControlPointMatrix.rows(); ++zeroOrderIdx) {
    controlPoint +=
      prefactor(idx, zeroOrderIdx, k, zerothOrderControlPointMatrix, zerothOrderKnotVector) * zerothOrderControlPointMatrix.row(zeroOrderIdx);
  }
  return controlPoint;
}


/* calculate matrices and vectors for the linear system of equations */
Eigen::MatrixXd BSplineMerger::calculateKV() {
  Eigen::MatrixXd KV(p_,p_);
  KV.setZero();
  double sum;

  for (unsigned k = 0; k <= p_-1; ++k) {
    for (unsigned i = m_-p_+1; i <= m_; ++i) {
      sum = 0;
      for (unsigned a = m_-p_; a <= m_-k; ++a) {
        sum += +prefactor(a, i, k, S_, Vk_[0]) * bsBasis_.evaluate(a, p_ - k, m_-k, Vk_[k], Vk_[0](m_ + 1));
      }
      KV(k,i-(m_-p_+1)) = sum;
    }
  }
  return KV;
}
Eigen::MatrixXd BSplineMerger::calculateKW() {
  Eigen::MatrixXd KW(p_,p_);
  KW.setZero();
  double sum;

  for (unsigned k = 0; k <= p_-1; ++k) {
    for (unsigned j = 0; j <= p_-1; ++j) {
      sum = 0;
      for (unsigned b = 0; b <= p_-k; ++b) {
        sum += -prefactor(b, j, k, T_, Wk_[0]) * bsBasis_.evaluate(b, p_ - k, o_-k, Wk_[k], Wk_[0](p_));
      }
      KW(k,j) = sum;
    }
  }
  return KW;
}
Eigen::MatrixXd BSplineMerger::calculateIV() {
  Eigen::MatrixXd IV(p_,p_);
  IV.setZero();
  double sum;

  for (unsigned i = m_-p_+1; i <= m_; ++i) {
    for (unsigned k = 0; k <= p_-1; ++k) {
      sum = 0;
      for (unsigned a = m_-p_; a <= m_-k; ++a) {
        sum += +prefactor(a, i, k, S_, Vk_[0]) * bsBasis_.evaluate(a, p_ - k, m_-k, Vk_[k], Vk_[0](m_ + 1));
      }
      IV(i-(m_-p_+1),k) = 0.5*sum;
    }
  }
  return IV;
}
Eigen::MatrixXd BSplineMerger::calculateJW() {
  Eigen::MatrixXd JW(p_,p_);
  JW.setZero();
  double sum;

  for (unsigned j = 0; j <= p_-1; ++j) {
    for (unsigned k = 0; k <= p_-1; ++k) {
      sum = 0;
      for (unsigned b = 0; b <= p_-k; ++b) {
        sum += -prefactor(b, j, k, T_, Wk_[0]) * bsBasis_.evaluate(b, p_ - k, o_-k, Wk_[k], Wk_[0](p_));
      }
      JW(j, k) = 0.5*sum;
    }
  }
  return JW;
}
Eigen::MatrixXd BSplineMerger::calculateGV(){
  Eigen::MatrixXd GV(mg_+1,p_);
  GV.setZero();

  for (unsigned g = 0; g <= mg_; ++g) {
    for (unsigned i = m_-p_+1; i <= m_; ++i) {
      GV(g,i-(m_-p_+1)) = +bsBasis_.evaluate(i, p_, m_, Vk_[0], constrainedParamsA_[g]);
    }
  }
  return GV;
}
Eigen::MatrixXd BSplineMerger::calculateHW(){
  Eigen::MatrixXd HW(oh_+1,p_);
  HW.setZero();

  for (unsigned h = 0; h <= oh_; ++h) {
    for (unsigned i = 0; i <= p_-1; ++i) {
      HW(h,i) = +bsBasis_.evaluate(i, p_, o_, Wk_[0], constrainedParamsB_[h]);
    }
  }
  return HW;
}
Eigen::MatrixXd BSplineMerger::calculateIpV() {
  Eigen::MatrixXd IpV(p_,mg_+1);
  IpV.setZero();

  for (unsigned i = m_-p_+1; i <= m_; ++i) {
    for (unsigned g = 0; g <= mg_; ++g) {
      IpV(i-(m_-p_+1), g) = -0.5* bsBasis_.evaluate(i, p_, m_, Vk_[0], constrainedParamsA_[g]);
    }
  }
  IpV *= 0.5;
  return  IpV;
}
Eigen::MatrixXd BSplineMerger::calculateJppW() {
  Eigen::MatrixXd JppW(p_,oh_+1);
  JppW.setZero();

  for (unsigned j = 0; j <= p_-1; ++j) {
    for (unsigned h = 0; h <= oh_; ++h) {
      JppW(j, h) = -0.5* bsBasis_.evaluate(j, p_, o_, Wk_[0], constrainedParamsB_[h]);
    }
  }
  JppW *= 0.5;
  return JppW;
}
Eigen::MatrixXd BSplineMerger::constructNmat() {
  Eigen::MatrixXd Nmat(3*p_+numberOfConstraintsA_+numberOfConstraintsB_,
                       3*p_+numberOfConstraintsA_+numberOfConstraintsB_);

  /* calculate block matrices for the p recise merge */
  Nmat.setZero();
  Nmat.topLeftCorner(2*p_,2*p_) = Eigen::MatrixXd::Identity(2*p_,2*p_);
  Nmat.block(2*p_,0 ,p_,p_) = calculateKV();;
  Nmat.block(2*p_,p_,p_,p_) = calculateKW();

  Nmat.block(0 ,2*p_,p_,p_) = calculateIV();
  Nmat.block(p_,2*p_,p_,p_) = calculateJW();

  /* calculate block matrices for the constraints */
  if (numberOfConstraintsA_>0) {
    mg_ = (unsigned) constrainedParamsA_.size() - 1;
    Nmat.block(3*p_,0 ,numberOfConstraintsA_,p_)= calculateGV();
    Nmat.block(0 ,3*p_,p_,numberOfConstraintsA_) = calculateIpV();
  }
  if (numberOfConstraintsB_>0) {
    oh_ = (unsigned) constrainedParamsB_.size() - 1;
    Nmat.block(3*p_+numberOfConstraintsA_,p_,numberOfConstraintsB_,p_)= calculateHW();
    Nmat.block(0 ,3*p_+numberOfConstraintsA_,p_,numberOfConstraintsB_) = calculateJppW();
  }

  //std::cout << Nmat << std::endl<< std::endl;

  return Nmat;
}
Eigen::MatrixXd BSplineMerger::calculateKconst() {
  Eigen::MatrixXd Kconst(p_,dim_);
  Kconst.setZero();

  Eigen::VectorXd sumA(dim_);
  Eigen::VectorXd sumB(dim_);

  for (unsigned k = 0; k <= p_-1; ++k) {
    sumA.setZero();
    sumB.setZero();
    for (unsigned i = m_-p_; i <= m_; ++i) {
      for (unsigned a = m_-p_; a <= m_-k; ++a) {
        sumA += prefactor(a, i, k, S_, Vk_[0]) * bsBasis_.evaluate(a, p_ - k, m_-k, Vk_[k], Vk_[0](m_ + 1)) * S_.row(i);
      }
    }
    for (unsigned j = 0; j <= p_; ++j) {
      for (unsigned b = 0; b <= p_-k; ++b) {
        sumB += -prefactor(b, j, k, T_, Wk_[0]) * bsBasis_.evaluate(b, p_ - k, o_-k, Wk_[k], Wk_[0](p_)) * T_.row(j);
      }
    }
    Kconst.row(k) = -sumA-sumB;
  }
  return Kconst;
}
Eigen::MatrixXd BSplineMerger::constructConstMat() {
  Eigen::MatrixXd Const (3*p_+numberOfConstraintsA_+numberOfConstraintsB_,dim_);
  Const.setZero();

  Const.block(2*p_,0,p_,dim_) = calculateKconst();

  return Const;
}
void BSplineMerger::initializeSolver(){
  //svdOfNmat_.compute(constructNmat());
  svdOfNmat_.compute(constructNmat(), Eigen::ComputeThinU | Eigen::ComputeThinV);
}
void BSplineMerger::solveLinearEquationSystem(){
  shifts_.resize(3*p_+numberOfConstraintsA_+numberOfConstraintsB_,dim_);
  shifts_.setZero();
  shifts_ = svdOfNmat_.solve(constructConstMat());
}
void BSplineMerger::generateShiftedControlPointVectorsOfSpline1and2(){
  Sshifted_ = S_;
  Sshifted_.bottomRows(p_) += shifts_.topRows(p_);

  Tshifted_ = T_;
  Tshifted_.topRows(p_) += shifts_.block(p_,0,p_,dim_);
}

Eigen::VectorXd BSplineMerger::reverseKnotVector(const Eigen::VectorXd U){
  Eigen::VectorXd Urev = (U.reverse().array()-1)*-1;
  return Urev;
}

Eigen::VectorXd BSplineMerger::adjustKnotVector(const Eigen::VectorXd & Uleft,
                                                const unsigned nLeft,
                                                const Eigen::VectorXd & Uright) {
  Eigen::VectorXd UleftAdj;
  UleftAdj.resize( (nLeft+2)+(p_) );
  UleftAdj.head(nLeft+2) = Uleft.head(nLeft+2);
  UleftAdj.tail(p_) = Uright.segment(p_+1,p_).array()+1;
  return UleftAdj;
}

void BSplineMerger::adjustKnotVectorsOfBothSplines(){
  Vadj_= adjustKnotVector(Vk_[0],m_,Wk_[0]);
  Vrev_ = reverseKnotVector(Vk_[0]);
  Wrev_ = reverseKnotVector(Wk_[0]);
  Wadjrev_ = adjustKnotVector(Wrev_,o_,Vrev_);
  Wadj_ = reverseKnotVector(Wadjrev_).array()+1;
}

void BSplineMerger::createMergedKnotVector(){
  /* construct the knot vector of the merged spline */
  Umerged_.resize(m_+2 + o_+1);
  Umerged_.head(m_+2) = Vadj_.head(m_+2);
  Umerged_.tail(o_+1) = Wadj_.tail(o_+1);

  /* rescale the knot vector of the merged spline to [0,1] by dividing through the last element */
  Umerged_ /= Umerged_(m_+o_+2);
}

void BSplineMerger::adjustParametrization(){
  /* TODO missing - can be necessary for arc-length parametrized splines (if the spline itself is arc-length paramtrized)
   * if the spline is arc-length parametrized by a second spline, then there is no need for adjusting the
   * parametrization at this stage */
}

Eigen::MatrixXd BSplineMerger::adjustShiftedControlPoints(const Eigen::MatrixXd & shiftedControlPoints,
                                                          const Eigen::VectorXd & originalKnotVector,
                                                          const Eigen::VectorXd & adjustedKnotVector,
                                                          const unsigned n){
  Eigen::MatrixXd Padj;
  Padj.resize(n+1,dim_);
  std::vector<std::vector<Eigen::VectorXd>> Qki(n+1);

  for (unsigned i = 0; i <= n-p_; ++i) {
    Qki[i].push_back(shiftedControlPoints.row(i));
  }
  for (unsigned k = 0; k <= p_-1; ++k) {
    Qki[n-p_+1].push_back( generateDerivativeControlPoint(n-p_+1,k,shiftedControlPoints,originalKnotVector) );
  }
  for (unsigned i = n-p_+2; i <= n; ++i) {
    Qki[i].resize(n-i+1);
    for (int k = n-i; k >=0; --k) {
      Qki[i][k] = ((adjustedKnotVector[i+p_]-adjustedKnotVector[i+k])/double(p_-k)) * Qki[i-1][k+1] + Qki[i-1][k];
    }
  }
  for (unsigned i = 0; i <= n; ++i) {
    Padj.row(i) = Qki[i][0];
  }
  return Padj;
}

void BSplineMerger::adjustShiftedControlPointsOfBothSplines() {

  /* obtain adjusted, shifted control points of the second spline */
  Sadj_ = adjustShiftedControlPoints(Sshifted_,Vk_[0],Vadj_,m_);

  /* reverse control points of the second spline */
  Eigen::MatrixXd P2shiftedrev = Tshifted_.colwise().reverse();

  /* obtain adjusted, shifted control points of the reversed second spline */
  Tadjrev_ = adjustShiftedControlPoints(P2shiftedrev,Wrev_,Wadjrev_,o_);

  /* obtain adjusted, shifted control points of the second spline */
  Tadj_ = Tadjrev_.colwise().reverse();
}

void BSplineMerger::generateControlPointVectorOfMergedSpline(){
  Pmerged_.resize(m_+1 + o_+1-p_,dim_);
  Pmerged_.topRows(m_+1) = Sadj_;
  Pmerged_.bottomRows(o_+2-p_) = Tadj_.bottomRows(o_+2-p_);
}

void BSplineMerger::printResultsForMathematica(){

  std::cout << "{Um,U1,U1adj,U2,U2adj,"
    "P1,P1shifted,P1adj,"
    "P2,P2shifted,P2adj,"
    "Pmerged}={";
  ContainerConverter::printVectorXdForMathematica(Umerged_);
  std::cout << "," << std::endl;
  ContainerConverter::printVectorXdForMathematica(Vk_[0]);
  std::cout << "," << std::endl;
  ContainerConverter::printVectorXdForMathematica(Vadj_);
  std::cout << "," << std::endl;
  ContainerConverter::printVectorXdForMathematica(Wk_[0]);
  std::cout << "," << std::endl;
  ContainerConverter::printVectorXdForMathematica(Wadj_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(S_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(Sshifted_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(Sadj_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(T_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(Tshifted_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(Tadj_);
  std::cout << "," << std::endl;
  ContainerConverter::printMatrixXdForMathematica(Pmerged_);
  std::cout << "};" << std::endl;
}