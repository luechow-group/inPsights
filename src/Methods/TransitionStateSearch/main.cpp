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
