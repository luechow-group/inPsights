//
// Created by Michael Heuer on 10.05.16.
//

#include "BSplineFromControlPolygon.h"

BSplineFromControlPolygon::BSplineFromControlPolygon(const Eigen::MatrixXd & dataPoints,
                                                     const unsigned degree,
                                                     const bool uniformKnotVector)
  : BSplineGenerator(dataPoints, degree), uniformKnotVector_(uniformKnotVector)
{
  initializeGenerator();
}

void BSplineFromControlPolygon::initializeGenerator(){
  if(uniformKnotVector_) {
    generateKnotVectorByUniformMethod();
  }
  else {
    generateParametersByChordLengthMethod();
    generateKnotVectorByDeBoorsMethod();
  }
  generateControlPointMatrix();
}

void BSplineFromControlPolygon::generateControlPointMatrix() {
  P_ = R_;
}

