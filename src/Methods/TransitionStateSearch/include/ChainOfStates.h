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

#ifndef INPSIGHTS_CHAINOFSTATES_H
#define INPSIGHTS_CHAINOFSTATES_H

#include <Eigen/Core>

class ChainOfStates{
public:
    ChainOfStates(const Eigen::MatrixXd& coordinates, const Eigen::VectorXd& values)
            : coordinates_(coordinates),
              values_(values)
    {}
    ChainOfStates(const Eigen::MatrixXd& coordinates)
            : coordinates_(coordinates),
              values_(this->statesNumber())
    {}

    long statesNumber(){return coordinates_.cols();}
    long coordinatesNumber(){return coordinates_.rows();}

    const Eigen::MatrixXd & coordinates() { return coordinates_; }
    const Eigen::VectorXd & values() { return values_; };
    //Eigen::MatrixXd   coordinatesCopy() { return coordinates_; }
    //const Eigen::MatrixXd & coordinates() { return coordinates_; }

    Eigen::VectorXd coordinatesAsVector();

    void storeCoordinateVectorInChain(long coordinatesNumber, long statesNumber, Eigen::VectorXd &coordinateVector);

    void setCoordinates(Eigen::MatrixXd coordinates){
      coordinates_ = coordinates;
    }

    void setValues(const Eigen::VectorXd & values){
      values_ = values;
    }

private:
    //Eigen::MatrixXd coordinates_;
    Eigen::MatrixXd coordinates_;//TODO should be row major
    Eigen::VectorXd values_;
};


#endif //INPSIGHTS_CHAINOFSTATES_H
