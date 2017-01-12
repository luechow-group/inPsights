#ifndef QUATERNIONFIT_HPP
#define QUATERNIONFIT_HPP

//#define EIGEN_USE_MKL_ALL

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

/* doi: 10.1021/ct501155k */

class QuaternionFit {
public:
  explicit QuaternionFit(const Eigen::MatrixXd &refMat);

  QuaternionFit(const Eigen::MatrixXd &refMat,
                const Eigen::MatrixXd &fitMat);

  QuaternionFit(const Eigen::MatrixXd &refMat,
                const Eigen::MatrixXd &fitMat,
                const Eigen::VectorXd &weights);

  virtual ~QuaternionFit();

  void align();

  void align(const Eigen::MatrixXd &fitMat);

  void align(const Eigen::MatrixXd &fitMat, const Eigen::VectorXd &weights);

  void calcRotMat();

  Eigen::Matrix4d getCMat() const { return c_; };

  Eigen::Vector4d getEigenvalues() const { return eigenvalues_; };

  Eigen::Matrix4d getEigenvectors() const { return eigenvectors_; };

  double getMaxEigenvalue() const { return maxEigenvalue_; };

  Eigen::Quaternion<double> getQuaternion() const { return quaternion_; };

  Eigen::Matrix3d getRotationMatrix() const { return rotMat_; };

  Eigen::Vector3d getTransVector() const { return transVector_; };

  Eigen::Vector3d getRefCenter() const { return refCenter_; };

  Eigen::Vector3d getFitCenter() const { return fitCenter_; };

  Eigen::MatrixXd getAlignedMat() const { return fittedMat_; };

  Eigen::Vector3d calcCenter(const Eigen::MatrixXd &mat, const Eigen::VectorXd &weights);

  double getRotRMSD() const { return rotRMSD_; };

  double getRMSD() const;


private:

  Eigen::VectorXd weights_;

  Eigen::MatrixXd refMat_;
  Eigen::MatrixXd fitMat_;

  Eigen::Vector3d refCenter_;
  Eigen::Vector3d fitCenter_;
  Eigen::Vector3d transVector_;

  Eigen::Matrix4d c_;
  Eigen::Vector4d eigenvalues_;
  Eigen::Matrix4d eigenvectors_;
  double maxEigenvalue_;

  Eigen::Quaternion<double> quaternion_;
  Eigen::Matrix3d rotMat_;

  Eigen::MatrixXd fittedMat_;
  double rotRMSD_;
};

Eigen::Matrix3d rotMatFromQuat(const Eigen::Vector4d &q);

#endif  //QUATERNIONFIT_HPP
