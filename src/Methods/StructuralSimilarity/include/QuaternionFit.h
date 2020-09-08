// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_QUATERNIONFITNEW_H
#define INPSIGHTS_QUATERNIONFITNEW_H

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

class QuaternionFit {
public:
  explicit QuaternionFit(const Eigen::MatrixXd &referencePositionsVector);

  QuaternionFit(const Eigen::MatrixXd &refMat,
                const Eigen::MatrixXd &fitMat);

  QuaternionFit(const Eigen::MatrixXd &referencePositionsVector,
                const Eigen::MatrixXd &targetPositionsVector,
                const Eigen::VectorXd &weights);

  virtual ~QuaternionFit();
  
  Eigen::Vector3d calculateCenter(const Eigen::MatrixXd &mat, const Eigen::VectorXd &weights);

  Eigen::Matrix3d rotationMatrixFromQuaternion(const Eigen::Vector4d &q);
  
  void calculateCorrelationMatrix();
  void alignTargetWithReference();

  void align();

  void align(const Eigen::MatrixXd &targetPositionsVector);

  void align(const Eigen::MatrixXd &targetPositionsVector, const Eigen::VectorXd &weights);
  
  //double getRotRMSD() const { return rotRMSD_; };
  //double getRMSD() const;
  Eigen::MatrixXd getRotationMatrix(){ return rotationMatrix_; };
  

private:
  Eigen::VectorXd weights_;

  Eigen::MatrixXd referencePositionsVector_;
  Eigen::MatrixXd targetPositionsVector_;

  Eigen::Vector3d referenceCenter_;
  Eigen::Vector3d targetCenter_;
  Eigen::Vector3d translationVector_;

  Eigen::Matrix4d correlationMatrix_;
  Eigen::Vector4d eigenvalues_;
  Eigen::Matrix4d eigenvectors_;
  double maxEigenvalue_;

  Eigen::Quaternion<double> quaternion_;
  Eigen::Matrix3d rotationMatrix_;

  Eigen::MatrixXd fittedTargetPositionsVector_;
  double rotRMSD_;
  unsigned atomNumber_;

  void findMaximalEigenvalue();

  void calculateRMSD();//TODO fix this
};

#endif //INPSIGHTS_QUATERNIONFITNEW_H
