//
// Created by Michael Heuer on 10.05.16.
//

#include "BSpline.h"
#include "BSplineTools.h"

BSpline::BSpline(){}

BSpline::BSpline(const Eigen::VectorXd& knotVector,
                 const Eigen::MatrixXd& controlPoints,
                 const unsigned degree,
                 int precalculatedDerivatives) :
  p_(degree),
  dim_((unsigned) controlPoints.cols()),
  n_((unsigned) controlPoints.rows()-1),
  isReversed_(false)
{
  assert(p_ >= 1 && "A B-spline curve must have a degree of at least 1");
  assert(n_ >= p_ && "A B-spline curve must have more control points than its own degree.");
  assert( precalculatedDerivatives >= -1 );

  Uk_.clear();
  Pk_.clear();
  Uk_.push_back(knotVector);
  Pk_.push_back(controlPoints);

  if (precalculatedDerivatives == -1)
    highestCalculatedDerivative_ = p_;
  else if ( (0 <= precalculatedDerivatives) && (unsigned(precalculatedDerivatives) < p_) )
    highestCalculatedDerivative_ = unsigned(precalculatedDerivatives);

  calculateDerivatives(highestCalculatedDerivative_);
}

BSpline BSpline::getDerivativeBSpline(const unsigned k) {
  assert( k < p_ );
  return BSpline(deriveAndGetKnotVector(k), deriveAndGetControlPointMatrix(k),p_-k,0);
}

void BSpline::reverseSpline(){
  for (unsigned i = 0; i < Uk_.size(); ++i) {
    Eigen::VectorXd tempVec(Uk_[i]);
    Uk_[i]= (tempVec.reverse().array()-1)*(-1); }

  for (unsigned i = 0; i < Pk_.size(); ++i) {
    Eigen::MatrixXd tempMat(Pk_[i]);
    Pk_[i]= tempMat.colwise().reverse(); }

  isReversed_ = !isReversed_;
}

bool BSpline::isClampedAndNormed() const {
  Eigen::VectorXd head,tail;
  head = Eigen::VectorXd::Zero(p_+1);
  tail = Eigen::VectorXd::Ones(p_+1);

  return Uk_[0].head(p_+1).isApprox(head) && Uk_[0].tail(p_+1).isApprox(tail);
}

const Eigen::VectorXd& BSpline::getKnotVector(const unsigned k) const {
  assert( (k < p_) && "This derivative does not exist.");
  assert( (k <= highestCalculatedDerivative_) && "This derviative is not calculated yet. Use deriveAndGetControlKnotVector(k) instead.");
  return Uk_[k];
}

const Eigen::VectorXd& BSpline::deriveAndGetKnotVector(const unsigned k){
  assert( (k < p_) && "This derivative does not exist.");
  if(k > highestCalculatedDerivative_) calculateDerivatives(k);
  return getKnotVector(k);
}

const std::vector<Eigen::VectorXd>& BSpline::getKnotVectors() const { return Uk_; }

const Eigen::MatrixXd& BSpline::getControlPointMatrix(const unsigned k) const {
  assert( (k < p_) && "This derivative does not exist.");
  assert( (k <= highestCalculatedDerivative_) && "This derviative is not calculated yet. Use deriveAndGetControlPointMatrix(k) instead.");
  return Pk_[k];
}

const Eigen::MatrixXd& BSpline::deriveAndGetControlPointMatrix(const unsigned k){
  assert( (k <= p_) && "This derivative does not exist.");
  if(k > highestCalculatedDerivative_) calculateDerivatives(k);
  return getControlPointMatrix(k);
}

const std::vector<Eigen::MatrixXd>& BSpline::getControlPointMatrices() const { return Pk_; }

void BSpline::deriveKnotVectors(const unsigned highestDerivativeToCalculate) {
  assert( (highestDerivativeToCalculate >= 0) && (highestDerivativeToCalculate <= p_) );

  for (unsigned k = (unsigned)Uk_.size(); k <= highestDerivativeToCalculate; ++k) {
    /* make the dth derivative knot vector from the previous d-1th derivative knot vector segment by dropping the
     * first and last segment */
    Eigen::VectorXd segmentOfPreviousOrderKnotVector(Uk_[k-1].segment(1,Uk_[k-1].size()-2));
    Uk_.push_back(segmentOfPreviousOrderKnotVector);
  }
}

void BSpline::deriveControlPointMatrices(const unsigned highestDerivativeToCalculate) {
  assert( (highestDerivativeToCalculate >= 0) && (highestDerivativeToCalculate <= p_) );

  for (unsigned k = (unsigned)Pk_.size(); k <= highestDerivativeToCalculate; ++k) {
    Eigen::MatrixXd Pnew(n_-k+1,dim_);
    for (unsigned i = 0; i <= n_-k; ++i) {
      Pnew.row(i) = deriveControlPoint(i, k);
    }
    Pk_.push_back(Pnew);
  }
}

void BSpline::calculateDerivatives(const unsigned highestDerivativeToCalculate) {
  assert( highestDerivativeToCalculate >= 0 );

  if ( highestDerivativeToCalculate > p_) {
    deriveKnotVectors(highestDerivativeToCalculate);
    deriveControlPointMatrices(highestDerivativeToCalculate);
    highestCalculatedDerivative_ = p_;
  }
  else {
    deriveKnotVectors(highestDerivativeToCalculate);
    deriveControlPointMatrices(highestDerivativeToCalculate);
    highestCalculatedDerivative_ = highestDerivativeToCalculate;
  }
}

Eigen::VectorXd BSpline::deriveControlPoint(const unsigned i, const unsigned k) const {
  if (k == 0)
    return Pk_[0].row(i);
  else if ( Uk_[0](i+p_+1) == Uk_[0](i+k) )
    return Eigen::VectorXd(Eigen::VectorXd::Zero(dim_));
  else
    return double(p_-k+1) / (Uk_[0](i+p_+1) - Uk_[0](i +k)) *
           (deriveControlPoint(i+1,k-1) - deriveControlPoint(i,k-1));
}

Eigen::VectorXd BSpline::deriveAndEvaluateNaive(const double u, const unsigned k) {
  assert( k < p_ && "Derivative order has to be less than the B-Spline degree.");
  assert( (0.0 <= u) && (u <= 1.0) );

  if ( k > highestCalculatedDerivative_ ) calculateDerivatives(k);

  Eigen::VectorXd C(dim_);
  C.setZero();
  for (unsigned i = 0; i <= n_-k ; ++i) {
    C += bsBasis_.evaluate(i, p_-k, n_-k, Uk_[k], u)*Pk_[k].row(i);
  }
  return C;
}

Eigen::VectorXd BSpline::operator()(const double u, const unsigned k) const {
  return evaluate(u, k);
}

Eigen::VectorXd BSpline::deriveAndEvaluate(const double u, const unsigned k){
  assert( k >= 0 && "Derivative order has to be greater than or equal to zero.");
  assert( (0.0 <= u) && (u <= 1.0) );

  if (highestCalculatedDerivative_ < k) {
    if (k <= p_) calculateDerivatives(k);
    else calculateDerivatives(p_);
  }

  return evaluate(u, k);
}

Eigen::VectorXd BSpline::evaluate(const double u, const unsigned k) const {
  assert( k >= 0 && "Derivative order has to be greater than or equal to zero.");
  assert( (0.0 <= u) && (u <= 1.0) );

  if ( k > p_) return Eigen::VectorXd(Eigen::VectorXd::Zero(dim_));
  else {
    unsigned l = findIdxOfLowerOrEqualDomainKnot(u, k);
    return deBoorAlgorithm(u,l,p_-k,k);
  }
}

Eigen::VectorXd BSpline::deBoorAlgorithm(const double u,const unsigned i,const unsigned p,const unsigned k) const {
  Eigen::VectorXd result;
  if (p == 0) result = Pk_[k].row(i);
  else {
    double alpha = (u-Uk_[k](i))/(Uk_[k](i+(p_-k)+1-p)-Uk_[k](i));
    result = (1-alpha) * deBoorAlgorithm(u,i-1,p-1,k) + alpha * deBoorAlgorithm(u,i,p-1,k);
  }
  return result;
}

unsigned BSpline::findIdxOfLowerOrEqualDomainKnot(const double u, const unsigned k) const {
  assert( k <= p_);

  // start at the beginning of the domain [u_{p-k}^{(k)},u_{n+1-k}^{(k)}]
  unsigned i = p_-k;
  while ( (u >= Uk_[k](i+1)) && ( (i+1)<(n_+1-k)) ) ++i;
  return i;
}
