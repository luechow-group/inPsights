/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
