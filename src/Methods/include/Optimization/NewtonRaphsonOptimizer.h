//
// Created by Michael Heuer on 26.01.17.
//

#ifndef AMOLQCGUI_NEWTONRAPHSONMINIMIZATION_H
#define AMOLQCGUI_NEWTONRAPHSONMINIMIZATION_H

#include "NewtonTypeOptimizer.h"
//#include "HessianFromGradientCalculator.h"

//#include <cppoptlib/solver/newtondescentsolver.h>

class NewtonRaphsonOptimizer : NewtonTypeOptimizer {

public:
  NewtonRaphsonOptimizer();

  void performStep() override;
  void constructHessian() override;
  //TODO make it more general => just getHessian?
  // maybe switch between provide hessian and self-calculated one by enum
  // => wrapping alglib

private:
  Eigen::MatrixXd hessian_;
  //HessianFromGradientCalculator hessianFromGradientCalculator_;
};

#endif //AMOLQCGUI_NEWTONRAPHSONMINIMIZATION_H
