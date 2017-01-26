//
// Created by Michael Heuer on 23.01.17.
//

#ifndef AMOLQCGUI_OPTIMIZERFACTORY_H
#define AMOLQCGUI_OPTIMIZERFACTORY_H

#include <memory>
#include <vector>

#include "LBFGSOptimizer.h"

enum class OptimizationMethod_t {
  LBFGS
};

class OptimizerFactory{
public:
  static std::unique_ptr<IOptimizer> createOptimizer(OptimizationMethod_t optimizationMethod)
  {
    std::unique_ptr<IOptimizer> method;

    switch(optimizationMethod){
      case OptimizationMethod_t::LBFGS:
        method = std::make_unique<LBFGSOptimizer>();
        break;
      default:
        break;
    }
    return method;
  }
};

#endif //AMOLQCGUI_OPTIMIZERFACTORY_H
