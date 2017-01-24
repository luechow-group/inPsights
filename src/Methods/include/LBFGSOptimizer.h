//
// Created by Michael Heuer on 23.01.17.
//

#ifndef AMOLQCGUI_LBFGSOPTIMIZER_H
#define AMOLQCGUI_LBFGSOPTIMIZER_H

#include "IOptimizer.h"
//#include "ElectronicWaveFunction.h"
//#include <optimization.h> //alglib


class LBFGSOptimizer : public IOptimizer {

public:
  LBFGSOptimizer();

  void performStep() override;

  // LBFGS specific methods
  //std::vector<std::vector<CartesionCoordinate>> getHessian(){ return hessian_; };

private:
  /*void functionAndGradient(const alglib::real_1d_array &x,
                           double &func,
                           alglib::real_1d_array &grad, void *ptr);*/

  //alglib::minlbfgsstate state_;
  //alglib::minlbfgsreport rep_;
  //alglib::real_1d_array AX_;
  //ElectronicWaveFunction wf_;
  //std::vector<std::vector<CartesionCoordinate>> hessian_;
};

#endif //AMOLQCGUI_LBFGSOPTIMIZER_H
