//
// Created by Michael Heuer on 16.01.17.
//
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

#ifndef AMOLQCGUI_QUATERNIONFITNEW_H
#define AMOLQCGUI_QUATERNIONFITNEW_H

class QuaternionFit {
public:
  explicit QuaternionFit(const Eigen::MatrixXd &referencePositionCollection);

  QuaternionFit(const Eigen::MatrixXd &refMat,
                const Eigen::MatrixXd &fitMat);

  QuaternionFit(const Eigen::MatrixXd &referencePositionCollection,
                const Eigen::MatrixXd &targetPositionCollection,
                const Eigen::VectorXd &weights);

  virtual ~QuaternionFit();
  
  Eigen::Vector3d calculateCenter(const Eigen::MatrixXd &mat, const Eigen::VectorXd &weights);

  Eigen::Matrix3d rotationMatrixFromQuaternion(const Eigen::Vector4d &q);
  
  void calculateCorrelationMatrix();
  void alignTargetWithReference();

  void align();

  void align(const Eigen::MatrixXd &targetPositionCollection);

  void align(const Eigen::MatrixXd &targetPositionCollection, const Eigen::VectorXd &weights);
  
  //double getRotRMSD() const { return rotRMSD_; };
  //double getRMSD() const;
  Eigen::MatrixXd getRotationMatrix(){ return rotationMatrix_; };

  void checkSpecialCases();

private:
  Eigen::VectorXd weights_;

  Eigen::MatrixXd referencePositionCollection_;
  Eigen::MatrixXd targetPositionCollection_;

  Eigen::Vector3d referenceCenter_;
  Eigen::Vector3d targetCenter_;
  Eigen::Vector3d translationVector_;

  Eigen::Matrix4d correlationMatrix_;
  Eigen::Vector4d eigenvalues_;
  Eigen::Matrix4d eigenvectors_;
  double maxEigenvalue_;

  Eigen::Quaternion<double> quaternion_;
  Eigen::Matrix3d rotationMatrix_;

  Eigen::MatrixXd fittedTargetPositionCollection_;
  double rotRMSD_;
  unsigned atomNumber_;

  void findMaximalEigenvalue();

  void calculateRMSD();//TODO fix this
};

#endif //AMOLQCGUI_QUATERNIONFITNEW_H
