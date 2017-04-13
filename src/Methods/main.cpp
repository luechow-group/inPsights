//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include "ElectronicWaveFunctionProblem.h"
#include "StringMethodProblem.h"
#include "solver/bfgsnssolver.h"

int main(int argc, char const *argv[]) {

    //StringMethod<cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem>> stringMethod(4); // is a ProblemObserver

    Eigen::VectorXd xA(2*3);
    Eigen::VectorXd x1(2*3);
    Eigen::VectorXd x2(2*3);
    Eigen::VectorXd x3(2*3);
    Eigen::VectorXd xB(2*3);

    xA << 0.0,0.0,+0.7,0.0,0.0,-0.7;
    x1 << 0.0,0.0,+0.3,0.0,0.0,-0.3;
    x2 << 0.0,0.0,+0.0,0.0,0.0,+0.0;
    x3 << 0.0,0.0,-0.3,0.0,0.0,+0.3;
    xB << 0.0,0.0,-0.7,0.0,0.0,+0.7;

    Eigen::VectorXd systemVector(5*2*3);
    systemVector << xA, x1, x2, x3, xB;
    std::cout << systemVector;

    StringMethodProblem stringMethodProblem(5,2*3);
    std::cout << stringMethodProblem.value(systemVector) << std::endl;
}
