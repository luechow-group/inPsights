//
// Created by Michael Heuer on 26.01.17.
//

#ifndef AMOLQCGUI_NEWTONTYPEOPTIMIZER_H
#define AMOLQCGUI_NEWTONTYPEOPTIMIZER_H

#include <Eigen/Core>
#include "IOptimizer.h"

class NewtonTypeOptimizer : public IOptimizer {
public:
  NewtonTypeOptimizer()
    : IOptimizer() {

  }

  virtual void constructHessian() = 0;

private:
  Eigen::MatrixXd hessian_;
};

#endif //AMOLQCGUI_NEWTONTYPEOPTIMIZER_H
