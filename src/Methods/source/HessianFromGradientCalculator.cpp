//
// Created by Michael Heuer on 26.01.17.
//

#include "assert.h" //TODO why needed?
#include "HessianFromGradientCalculator.h"

HessianFromGradientCalculator::HessianFromGradientCalculator() {

}

void HessianFromGradientCalculator::constructHessianFromGradient() {

  // access to gradient calculating function is required since multiple gradient around the point of interest are
  // neeede to compute the hessian numberical by the finite difference method
  // => gradient ist the function itself => nested call
  // make numerical gradient calculator class and  access alglib for that

  assert(false && "to be implemented");
}
