//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include "ElectronicWaveFunctionProblem.h"
#include "StringOptimizationProblem.h"
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

    Eigen::MatrixXd initialChain(5,2*3);

    initialChain.row(0) = xA;
    initialChain.row(1) = x1;
    initialChain.row(2) = x2;
    initialChain.row(3) = x3;
    initialChain.row(4) = xB;

    StringMethod stringMethod(initialChain);

    std::cout << stringMethod.getChain().coordinates() << std::endl;

    stringMethod.optimizeString();


    /*
    stringMethod.reparametrizeString();
    stringMethod.discretizeSplineToChain();
    stringMethod.performStep();*/
}
