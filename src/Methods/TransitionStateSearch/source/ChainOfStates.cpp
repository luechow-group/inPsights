//
// Created by heuer on 12.05.17.
//

#include "ChainOfStates.h"

Eigen::VectorXd ChainOfStates::coordinatesAsVector() {
  // map is performed in column major
  return Eigen::Map<Eigen::VectorXd>(coordinates_.data(), statesNumber()*coordinatesNumber());
};


void ChainOfStates::storeCoordinateVectorInChain(long coordinatesNumber, long statesNumber,
                                                 Eigen::VectorXd &coordinateVector) {
  // map is performed in column major
  coordinates_ = Eigen::Map<Eigen::MatrixXd>(coordinateVector.data(), coordinatesNumber, statesNumber);
}
