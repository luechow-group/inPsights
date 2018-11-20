//
// Created by heuer on 12.05.17.
//

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
