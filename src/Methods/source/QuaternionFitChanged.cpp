#include "QuaternionFit.h"

#include <Eigen/Eigenvalues>

#include "iostream"
using namespace std;
using namespace Eigen;

QuaternionFit::QuaternionFit(const MatrixXd& refMat,
                             const MatrixXd& fitMat,
                             const VectorXd& weights)
	: weights_(weights),
		refMat_(refMat),
		fitMat_(fitMat)
{
	assert(refMat_.rows() == fitMat_.rows());
	assert(refMat_.rows() == weights_.rows());

	align();
}

QuaternionFit::QuaternionFit(const MatrixXd& refMat,
														 const MatrixXd& fitMat)
	: refMat_(refMat),
		fitMat_(fitMat)
{
	assert(refMat.rows() == fitMat.rows());
	weights_ = VectorXd::Ones(refMat_.rows());

	align();
}

QuaternionFit::QuaternionFit(const MatrixXd& refMat)
	: refMat_(refMat)
 {
	weights_ = VectorXd::Ones(refMat_.rows());
}

QuaternionFit::~QuaternionFit() {}

void QuaternionFit::align() {
	calcRotMat();
}

void QuaternionFit::align(const MatrixXd& fitMat) {
	assert(fitMat.cols() == 3);
	assert(fitMat.rows() == refMat_.rows());
	fitMat_ = fitMat;
	calcRotMat();
}

void QuaternionFit::align(const MatrixXd& fitMat, const VectorXd& weights) {
	assert(fitMat.cols() == 3);
	assert(fitMat.rows() == refMat_.rows());
	assert(refMat_.rows() == weights.rows());
	weights_ = weights;
	calcRotMat();
}

void QuaternionFit::calcRotMat() {
	int NAtoms = static_cast<int>(refMat_.rows());

	c_ = Matrix4d::Zero();

	refCenter_ = calcCenter(refMat_, weights_);
	fitCenter_ = calcCenter(fitMat_, weights_);
	transVector_ = fitCenter_ - refCenter_;

	//for (int i = 0; i < NAtoms; i++) {
	//	cRefMat.row(i) -= refCenter_;
	//	cFitMat.row(i) -= fitCenter_;
	//}

	double R11 = 0.0, R22 = 0.0, R33 = 0.0;
	double R23 = 0.0, R32 = 0.0;
	double R13 = 0.0, R31 = 0.0;
	double R12 = 0.0, R21 = 0.0;

	auto w = 0.0;

	for (int i = 0; i < NAtoms; i++) {
		auto fitPos = (fitMat_.row(i) - fitCenter_.transpose()).eval();
		auto refPos = (refMat_.row(i) - refCenter_.transpose()).eval();
		w = weights_(i);

		R11 += fitPos(0) * refPos(0) * w;
		R12 += fitPos(0) * refPos(1) * w;
		R13 += fitPos(0) * refPos(2) * w;
		R21 += fitPos(1) * refPos(0) * w;
		R22 += fitPos(1) * refPos(1) * w;
		R23 += fitPos(1) * refPos(2) * w;
		R31 += fitPos(2) * refPos(0) * w;
		R32 += fitPos(2) * refPos(1) * w;
		R33 += fitPos(2) * refPos(2) * w;
	}

	c_(0, 0) = R11 + R22 + R33;

	c_(0, 1) = R23 - R32;
	c_(1, 1) = R11 - R22 - R33;

	c_(0, 2) = R31 - R13;
	c_(1, 2) = R12 + R21;
	c_(2, 2) = -R11 + R22 - R33;

	c_(0, 3) = R12 - R21;
	c_(1, 3) = R31 + R13;
	c_(2, 3) = R23 + R32;
	c_(3, 3) = -R11 - R22 + R33;

	for (unsigned i = 0; i < 3; i++) {
		for (unsigned j = i+1; j < 4; j++) {
			c_(j, i) = c_(i, j);
		}
  }

  std::cout << std::endl << c_ << std::endl;

	SelfAdjointEigenSolver<Matrix4d> eigensolver(c_);
	eigenvalues_ = eigensolver.eigenvalues();
	eigenvectors_ = eigensolver.eigenvectors();

  // The eigenvector corresponding to the largest eigenvalue is the unit quaternion.
  // It can be used to construct the corresponding rotation matrix.
  int index;
	if (eigenvalues_.maxCoeff() > -eigenvalues_.minCoeff())
	{
		eigenvalues_.maxCoeff(&index);
		rotMat_ = rotMatFromQuat(eigenvectors_.col(index));
	} 
	else 
	{
		eigenvalues_.minCoeff(&index);
		rotMat_ = -rotMatFromQuat(eigenvectors_.col(index));
	}
	maxEigenvalue_ = eigenvalues_(index);

  std::cout << std::endl << rotMat_ << std::endl;

	fittedMat_ = fitMat_;
	for (unsigned i = 0; i < fittedMat_.rows(); i++) {
		fittedMat_.row(i) -= fitCenter_; 
		fittedMat_.row(i) = rotMat_ * fittedMat_.row(i).transpose();
		fittedMat_.row(i) += refCenter_;
	}

	auto sum = 0.0;
	for (int i = 0; i < NAtoms; i++) {
		auto fitPos = fitMat_.row(i) - fitCenter_.transpose();
		auto refPos = refMat_.row(i) - refCenter_.transpose();
		sum += fitPos.squaredNorm() + refPos.squaredNorm();
	}
	sum -= 2.0 * abs(maxEigenvalue_);

	if (sum < 1e-7) {
		rotRMSD_ = 0.0;
	} else {
		rotRMSD_ = sqrt(sum/NAtoms);
	}
}

Vector3d QuaternionFit::calcCenter(const MatrixXd& mat, const VectorXd& weights) {
	assert(weights.rows() == mat.rows());

	Vector3d center = Vector3d::Zero();
	for (unsigned i = 0; i < mat.rows(); i++)
		center += mat.row(i).transpose() * weights(i);

	center /= weights.sum();

	return center;
}

double QuaternionFit::getRMSD() const {
	auto sum = 0.0;

	for (unsigned i = 0; i < refMat_.rows(); i++)
		sum += (refMat_.row(i) - fittedMat_.row(i)).squaredNorm();

	if (sum < 1e-7) return 0.0;

	return sqrt(sum/(refMat_.rows()*3));
}

Matrix3d rotMatFromQuat(const Vector4d& q) {
	Matrix3d rotMat = Matrix3d::Zero();

	double q00, q11, q22, q33;
	double q01, q02, q03, q12, q13, q23;

	q00 = q(0) * q(0);
	q11 = q(1) * q(1);
	q22 = q(2) * q(2);
	q33 = q(3) * q(3);

	q01 = q(0) * q(1);
	q02 = q(0) * q(2);
	q03 = q(0) * q(3);
	q12 = q(1) * q(2);
	q13 = q(1) * q(3);
	q23 = q(2) * q(3);

  rotMat(0, 0) = (q00 + q11 - q22 - q33) * 0.5;
  rotMat(1, 0) = q12 - q03;
  rotMat(2, 0) = q13 + q02;

  rotMat(0, 1) = q12 - q03;
  rotMat(1, 1) = (q00 - q11 + q22 - q33) * 0.5;
  rotMat(2, 1) = q23 - q01;

  rotMat(0, 2) = q13 + q02;
  rotMat(1, 2) = q23 - q01;
  rotMat(2, 2) = (q00 - q11 - q22 + q33) * 0.5;

/*rotMat(0, 0) = q00 + q11 - q22 - q33;
	rotMat(1, 0) = 2.0 * (q12 + q03);
	rotMat(2, 0) = 2.0 * (q13 - q02);

	rotMat(0, 1) = 2.0 * (q12 - q03);
	rotMat(1, 1) = q00 - q11 + q22 - q33;
	rotMat(2, 1) = 2.0 * (q23 + q01);

	rotMat(0, 2) = 2.0 * (q13 + q02);
	rotMat(1, 2) = 2.0 * (q23 - q01);
	rotMat(2, 2) = q00 - q11 - q22 + q33;*/

	return rotMat;
}