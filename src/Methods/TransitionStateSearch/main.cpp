//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include "StringMethod.h"

int main(int argc, char const *argv[]) {

    //StringMethod<cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem>> stringMethod(4); // is a ProblemObserver

    Eigen::VectorXd xA(2*3);
    Eigen::VectorXd xB(2*3);

    xA << 0.0,  0.01,  0.70014273,  0.02, 0.09,-0.70014273;//0.0,0.00,+0.7,0.0,0.0,-0.7;
    xB << 0.01,-0.0, +0.70014273,  -0.05, 0.01,+0.70014273;//0.0,0.00,-0.7,0.0,0.0,+0.7;

  unsigned numberOfStates = 15;

  Eigen::MatrixXd initialChain(2*3,numberOfStates);

  Eigen::VectorXd delta = xB-xA;
  for (int i = 0; i < numberOfStates ; ++i) {
    double rel = double(i)/double(numberOfStates-1);

    initialChain.col(i) = xA + (delta * rel);
  }




    StringMethod stringMethod(initialChain);

    std::cout << stringMethod.getChain().coordinates() << std::endl;

    stringMethod.optimizeString();
}
