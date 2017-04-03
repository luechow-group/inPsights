//
// Created by Michael Heuer on 06.05.16.
//

#ifndef RTQC_CONTAINERCONVERTER_H
#define RTQC_CONTAINERCONVERTER_H

#include <Eigen/Core>
#include <fstream>
//#include <QList>

class BSpline;

/*! Contains static methods for the conversion of Delib types into Eigen vectors and matrices
 * mainly needed to convert PositionCollections and associated energies into a Eigen matrices to use them with
 * the B-Spline methods.
 * Also contains methods to print and write matrices, vectors, and B-splines.
 * */
class ContainerConverter {
public:
  using VectorNu = Eigen::Matrix<unsigned,Eigen::Dynamic,1>; // for element types
  using Vector3Nd = Eigen::VectorXd;

  // Matrix storing an exploration path of X points constituted by 3N atomic coordinates as doubles
  using MatrixX3Nd = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;

  // Matrix storing an exploration path of X points constituted by the energy E and 3N atomic coordinates as doubles
  using MatrixXE3Nd = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;

  ContainerConverter()= delete;


  /*
  static std::pair<Eigen::VectorXd,MatrixX3Nd> QListToVectorXdAndMatrixX3NdPair(const QList<PathElement> &qList);
  static QList<PathElement> QListFromVectorXdAndMatrixX3Nd(const Eigen::VectorXd & vectorXd,
                                                           const MatrixX3Nd & matrixX3Nd);

  static QList<PathElement> QListFromMatrixXE3Nd(const MatrixXE3Nd & qList);
  static MatrixXE3Nd QListToMatrixXE3Nd(const QList<PathElement> & qList);
  */

  static ContainerConverter::MatrixXE3Nd MatrixXE3NdFromVectorXdAndMatrixX3Nd(const Eigen::VectorXd & vectorXd,
                                                                              const MatrixX3Nd & matrixX3Nd);



  static void printVectorXdForMathematica(const Eigen::VectorXd &vec);
  static void printMatrixXdForMathematica(const Eigen::MatrixXd &mat);

  static void writeVectorXdForMathematica(const Eigen::VectorXd &vec, std::ofstream &fout);
  static void writeMatrixXdForMathematica(const Eigen::MatrixXd &mat, std::ofstream &fout);
};

#endif //RTQC_CONTAINERCONVERTER_H
