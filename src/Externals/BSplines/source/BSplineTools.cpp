//
// Created by Michael Heuer on 17.06.16.
//

#include "BSplineTools.h"

BSplineTools::BSplineTools(){
}

unsigned BSplineTools::findIdxOfLeftOrEqualDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U) {
  unsigned i = p;

  while ( (u >= U(i+1)) && ( i+1 < U.size()-p) ) {
    ++i;
    if( U(i+1)==U(i)) break; //stop iteration
  }
  return i;
}

unsigned BSplineTools::findIdxOfLeftDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U) {
  unsigned i = p;
  while ( (u > U(i+1)) && (i+1 < U.size()-p) ) ++i;
  return i;
}

unsigned BSplineTools::findIdxOfRightDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U) {
  unsigned i = (unsigned) U.size()-1-p;
  while ( (u < U(i-1)) && (i-1  > p ) ) --i;
  return i;
}

unsigned BSplineTools::findIdxOfRightOrEqualDomainKnot(const double u, const unsigned p, const Eigen::VectorXd &U) {
  unsigned i = (unsigned) U.size()-1-p;
  while ( (u <= U(i-1)) && (i-1  >= p ) ) {
    --i;
    if( U(i-1)==U(i)) break; //stop iteration
  }
  return i;
}

double BSplineTools::knotAverage(const unsigned i,
                                 const unsigned p,
                                 const Eigen::VectorXd &U){
  double knotSum = 0.0;
  for (unsigned j = i+1 ; j <= i+p; ++j) {
    knotSum += U(j);
  }
  return knotSum/double(p);
}

void BSplineTools::normalizeKnotVector(Eigen::VectorXd &U) {
  std::pair<double,double> oldLim = {U.minCoeff(),U.maxCoeff()};

  rescaleKnotVector(U,oldLim,{0,1});
}

Eigen::VectorXd BSplineTools::normalizedKnotVector(const Eigen::VectorXd &U) {
  std::pair<double,double> oldLim = {U.minCoeff(),U.maxCoeff()};

  return rescaledKnotVector(U,oldLim,{0,1});
}


double BSplineTools::rescaledKnot(const double knot,
                                 const std::pair<double,double> oldLim,
                                 const std::pair<double,double> newLim){

  double knotCopy = knot;
  knotCopy -= oldLim.first;
  knotCopy /= (oldLim.second-oldLim.first);
  knotCopy *= (newLim.second-newLim.first);
  knotCopy += newLim.first;
  return knotCopy;
}

Eigen::VectorXd BSplineTools::rescaledKnotVector(const Eigen::VectorXd &U,
                                                const std::pair<double,double> oldLim,
                                                const std::pair<double,double> newLim){
  auto Ucopy = U;
  rescaleKnotVector(Ucopy,oldLim,newLim);
  return Ucopy;
}


void BSplineTools::rescaleKnotVector(Eigen::VectorXd &U,
                                     const std::pair<double,double> oldLim,
                                     const std::pair<double,double> newLim){
  U.array() -= oldLim.first;
  U /= (oldLim.second-oldLim.first);
  U *= (newLim.second-newLim.first);
  U.array() += newLim.first;
}

//finite difference matrix for penalizing BSplines
int BSplineTools::differenceOperator(unsigned i, unsigned j, unsigned kappa) {
  if (kappa > 1){
    return  differenceOperator(i+1,j,kappa-1)-differenceOperator(i,j,kappa-1);
  }
  else if (kappa == 1){
    if( i+1 == j)  return 1;
    else if ( i == j) return -1;
    else return 0;
  }
  else return 0;
}

