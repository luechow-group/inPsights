//
// Created by heuer on 06.04.17.
//


#include "StringMethod.h"
#include "problem.h"

void StringMethod::resetString() {

}

void StringMethod::evaluateString(Eigen::VectorXd &x) {
    double f = problemReference.value(x);
    std::cout << f << std::endl;
}

void StringMethod::stepPerformed() {
   std::cout << problemReference.getObserverCount() << " oberservers are listening" << std::endl;
}
