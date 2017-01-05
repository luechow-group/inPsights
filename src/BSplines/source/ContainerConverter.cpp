//
// Created by Michael Heuer on 06.05.16.
//

#include <iostream>
#include "ContainerConverter.h"
#include "BSpline.h"

using namespace Eigen;
using namespace std;

/*
std::pair<Eigen::VectorXd,Eigen::MatrixXd>
ContainerConverter::QListToVectorXdAndMatrixX3NdPair(const QList<PathElement> &qList){

  unsigned nCoords = qList[0].second.size()*3;
  ContainerConverter::MatrixX3Nd matrixX3Nd(qList.size(),nCoords);
  Eigen::VectorXd vectorXd (qList.size());

  for (unsigned i = 0; i < qList.size(); ++i) {
    matrixX3Nd.row(i) = PositionCollectionToVector3Nd(qList[i].second);

    vectorXd(i) = qList[i].first;
  }
  return std::make_pair(vectorXd,matrixX3Nd);
}

QList<PathElement>
ContainerConverter::QListFromVectorXdAndMatrixX3Nd(const Eigen::VectorXd & vectorXd,
                                                   const MatrixX3Nd & matrixX3Nd){
  assert(matrixX3Nd.rows() == vectorXd.size());
  QList<PathElement> qList;

  for (unsigned i = 0; i < matrixX3Nd.rows(); ++i) {
    qList.push_back(std::make_pair(vectorXd(i),PositionCollectionFromVector3Nd(matrixX3Nd.row(i))));
  }

  return qList;
}

ContainerConverter::MatrixXE3Nd ContainerConverter::QListToMatrixXE3Nd(const QList<PathElement> &qList) {

  auto energyVectorAndCoordMatrixPair = QListToVectorXdAndMatrixX3NdPair(qList);

  unsigned nPoints = unsigned(energyVectorAndCoordMatrixPair.second.rows());
  unsigned nCoords = unsigned(energyVectorAndCoordMatrixPair.second.cols());

  MatrixXE3Nd energyCoordsMatrix(nPoints,1+nCoords);

  energyCoordsMatrix.block(0,0,nPoints,1) = energyVectorAndCoordMatrixPair.first;
  energyCoordsMatrix.block(0,1,nPoints,nCoords) = energyVectorAndCoordMatrixPair.second;

  return energyCoordsMatrix;
}

*/


ContainerConverter::MatrixXE3Nd
ContainerConverter::MatrixXE3NdFromVectorXdAndMatrixX3Nd(const Eigen::VectorXd & vectorXd,
                                                         const MatrixX3Nd & matrixX3Nd){
  assert(matrixX3Nd.rows() == vectorXd.size());

  ContainerConverter::MatrixXE3Nd matrixXE3Nd(matrixX3Nd.rows(),1+matrixX3Nd.cols());
  matrixXE3Nd.col(0) = vectorXd;
  matrixXE3Nd.rightCols(matrixX3Nd.cols()) = matrixX3Nd;

  return matrixXE3Nd;
}


void ContainerConverter::printMatrixXdForMathematica(const Eigen::MatrixXd &mat) {

  std::cout << "{";
  for (unsigned i = 0; i < mat.rows(); ++i) {
    std::cout << "{";
    for (unsigned j = 0; j < mat.cols(); ++j) {
      std::cout << mat(i,j);
      if(j < mat.cols()-1){ std::cout << ",";}
    }
    std::cout << "}";
    if(i < mat.rows()-1){ std::cout << ",";}
  }
  std::cout << "}";
}

void ContainerConverter::printVectorXdForMathematica(const Eigen::VectorXd &vec) {
  std::cout << "{";
  for (unsigned i = 0; i < vec.size(); ++i) {
      std::cout << vec(i);
    if(i < vec.rows()-1){ std::cout << ",";}
  }
  std::cout << "}";
}

void ContainerConverter::writeVectorXdForMathematica(const Eigen::VectorXd &vec, std::ofstream &fout) {
  fout << "";
  for (unsigned i = 0; i < vec.size(); ++i) {
    fout << vec(i);
    if(i < vec.rows()-1){ fout << ",";}
  }
  fout << "";
}

void ContainerConverter::writeMatrixXdForMathematica(const Eigen::MatrixXd &mat, std::ofstream &fout) {

  fout << "";
  for (unsigned i = 0; i < mat.rows(); ++i) {
    fout << "{";
    for (unsigned j = 0; j < mat.cols(); ++j) {
      fout << mat(i, j);
      if (j < mat.cols() - 1) { fout << ","; }
    }
    fout << "}";
    if (i < mat.rows() - 1) { fout << ","; }
  }
  fout << "";
}
