//
// Created by Michael Heuer on 16.01.17.
//

#include "QuaternionFit.h"
#include <Eigen/Eigenvalues>

#include "iostream"

using namespace std;
using namespace Eigen;

QuaternionFit::QuaternionFit(const MatrixXd& referencePositionsVector,
                             const MatrixXd& targetPositionsVector,
                             const VectorXd& weights)
  : weights_(weights),
    referencePositionsVector_(referencePositionsVector),
    targetPositionsVector_(targetPositionsVector),
    atomNumber_((unsigned)referencePositionsVector.rows())
{
  assert(atomNumber_ > 0);
  assert(referencePositionsVector_.rows() == targetPositionsVector_.rows());
  assert(referencePositionsVector_.rows() == weights_.rows());
  align();
}

QuaternionFit::QuaternionFit(const MatrixXd& referencePositionsVector,
                             const MatrixXd& targetPositionsVector)
  : referencePositionsVector_(referencePositionsVector),
    targetPositionsVector_(targetPositionsVector),
    atomNumber_((unsigned)referencePositionsVector.rows())
{
  assert(atomNumber_ > 0);
  assert(referencePositionsVector.rows() == targetPositionsVector.rows());

  weights_ = VectorXd::Ones(referencePositionsVector_.rows());

  align();
}

QuaternionFit::QuaternionFit(const MatrixXd& referencePositionsVector)
  : referencePositionsVector_(referencePositionsVector),
    atomNumber_((unsigned)referencePositionsVector.rows())
{
  assert(atomNumber_ > 0);
  weights_ = VectorXd::Ones(referencePositionsVector_.rows());
}

QuaternionFit::~QuaternionFit() {}

void QuaternionFit::align(const MatrixXd& targetPositionsVector, const VectorXd& weights) {
  assert(referencePositionsVector_.rows() == weights.rows());
  weights_ = weights;

  align(targetPositionsVector);
}

void QuaternionFit::align(const MatrixXd& targetPositionsVector) {
  assert(targetPositionsVector.cols() == 3);
  assert(targetPositionsVector.rows() == referencePositionsVector_.rows());
  targetPositionsVector_ = targetPositionsVector;

  align();
}

void QuaternionFit::align() {

  referenceCenter_ = calculateCenter(referencePositionsVector_, weights_);
  targetCenter_ = calculateCenter(targetPositionsVector_, weights_);

  assert(atomNumber_ > 0);

  if (atomNumber_ > 3) {
    calculateCorrelationMatrix();
    findMaximalEigenvalue();
    alignTargetWithReference();
    //calculateRMSD();
  }
  else if ( atomNumber_== 3){
    assert(false && "implementation needed");
  }
  else if ( atomNumber_== 2){
    // make target and reference collinear
    assert(false && "implementation needed");
  }
  else {
    // translate target to reference
    assert(false && "implementation needed");
  }
}

Eigen::Vector3d QuaternionFit::calculateCenter(const MatrixXd& mat, const VectorXd& weights) {
  assert(weights.rows() == mat.rows());

  Vector3d center = Vector3d::Zero();
  for (unsigned i = 0; i < mat.rows(); i++)
    center += mat.row(i).transpose() * weights(i);

  center /= weights.sum();

  return center;
}

void QuaternionFit::calculateCorrelationMatrix() {
  correlationMatrix_ = Matrix4d::Zero();
  double c11 = 0.0, c22 = 0.0, c33 = 0.0;
  double c23 = 0.0, c32 = 0.0;
  double c13 = 0.0, c31 = 0.0;
  double c12 = 0.0, c21 = 0.0;

  auto w = 0.0;

  for (unsigned i = 0; i < atomNumber_; i++) {

    // translate target and reference to origin
    auto targetPosition = (targetPositionsVector_.row(i) - targetCenter_.transpose()).eval();
    auto referencePosition = (referencePositionsVector_.row(i) - referenceCenter_.transpose()).eval();
    w = weights_(i);

    c11 += targetPosition(0) * referencePosition(0) * w;
    c12 += targetPosition(0) * referencePosition(1) * w;
    c13 += targetPosition(0) * referencePosition(2) * w;
    c21 += targetPosition(1) * referencePosition(0) * w;
    c22 += targetPosition(1) * referencePosition(1) * w;
    c23 += targetPosition(1) * referencePosition(2) * w;
    c31 += targetPosition(2) * referencePosition(0) * w;
    c32 += targetPosition(2) * referencePosition(1) * w;
    c33 += targetPosition(2) * referencePosition(2) * w;
  }

  correlationMatrix_ = Matrix4d::Zero();

  correlationMatrix_(0, 0) = c11 + c22 + c33;

  correlationMatrix_(0, 1) = c23 - c32;
  correlationMatrix_(1, 1) = c11 - c22 - c33;

  correlationMatrix_(0, 2) = c31 - c13;
  correlationMatrix_(1, 2) = c12 + c21;
  correlationMatrix_(2, 2) = -c11 + c22 - c33;

  correlationMatrix_(0, 3) = c12 - c21;
  correlationMatrix_(1, 3) = c13 + c31;
  correlationMatrix_(2, 3) = c23 + c32;
  correlationMatrix_(3, 3) = -c11 - c22 + c33;

  // mirror correlation matrix along diagonal
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = i+1; j < 4; j++) {
      correlationMatrix_(j, i) = correlationMatrix_(i, j);
    }
  }
}

void QuaternionFit::findMaximalEigenvalue() {
  SelfAdjointEigenSolver<Matrix4d> eigensolver(correlationMatrix_);
  eigenvalues_ = eigensolver.eigenvalues();
  eigenvectors_ = eigensolver.eigenvectors();


  int index;
  if (eigenvalues_.maxCoeff() >= -eigenvalues_.minCoeff()) {
    eigenvalues_.maxCoeff(&index);
    rotationMatrix_ = rotationMatrixFromQuaternion(eigenvectors_.col(index));
  } else {
    eigenvalues_.minCoeff(&index);
    rotationMatrix_ = -rotationMatrixFromQuaternion(eigenvectors_.col(index));
  }
  maxEigenvalue_ = eigenvalues_(index);
  //TODO check degeneracy
  //TODO select maximal, degenerate eigenvector with most zeros
}

Matrix3d QuaternionFit::rotationMatrixFromQuaternion(const Vector4d& q) {
  Matrix3d rotMat = Matrix3d::Zero();

  double q11, q22, q33;
  double q01, q02, q03, q12, q13, q23;

  q11 = q(1) * q(1);
  q22 = q(2) * q(2);
  q33 = q(3) * q(3);

  q01 = q(0) * q(1);
  q02 = q(0) * q(2);
  q03 = q(0) * q(3);
  q12 = q(1) * q(2);
  q13 = q(1) * q(3);
  q23 = q(2) * q(3);

  rotMat(0, 0) = 1 - 2*(q22 + q33);
  rotMat(1, 0) = 2*(q12 + q03);
  rotMat(2, 0) = 2*(q13 - q02);

  rotMat(0, 1) = 2*(q12 - q03);
  rotMat(1, 1) = 1 - 2*(q11 + q33);
  rotMat(2, 1) = 2*(q23 + q01);

  rotMat(0, 2) = 2*(q13 + q02);
  rotMat(1, 2) = 2*(q23 - q01);
  // PAPER IS WRONG! HERE IT IS CORRECT
  //http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/
  return rotMat;
}

void QuaternionFit::alignTargetWithReference() {

  fittedTargetPositionsVector_ = targetPositionsVector_;
  for (unsigned i = 0; i < fittedTargetPositionsVector_.rows(); i++) {
    // move to origin
    fittedTargetPositionsVector_.row(i) -= targetCenter_;
    // rotate
    fittedTargetPositionsVector_.row(i) = rotationMatrix_ * fittedTargetPositionsVector_.row(i).transpose();
    // move to reference center
    fittedTargetPositionsVector_.row(i) += referenceCenter_;
  }
}

void QuaternionFit::calculateRMSD() {
  auto sum = 0.0;

  // TODO FIX THIS

  for (unsigned i = 0; i < atomNumber_; i++) {
    auto targetPosition = targetPositionsVector_.row(i) - targetCenter_.transpose();
    auto referencePosition = referencePositionsVector_.row(i) - referenceCenter_.transpose();
    sum +=  targetPosition.squaredNorm() - referencePosition.squaredNorm(); // MACHT SO KEINEN SINN
  }
  sum -= 2.0 * abs(maxEigenvalue_);

  /*
  if (sum < 1e-7) {
    rotRMSD_ = 0.0;
  } else {
    rotRMSD_ = sqrt(sum / atomNumber_);
  }*/
  rotRMSD_ = sqrt(sum / atomNumber_);
}


