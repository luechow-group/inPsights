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
