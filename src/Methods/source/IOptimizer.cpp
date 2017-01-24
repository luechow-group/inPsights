//
// Created by Michael Heuer on 23.01.17.
//

#include <assert.h> //TODO why are several assert.h available here and in BSpline there is not?
#include "IOptimizer.h"

IOptimizer::IOptimizer()
  : stepCount_(0) {

}

void IOptimizer::iterate() {
    while( !converged() || (stepCount_<5) ) {
      stepCount_++;
      fetchGradient();
      performStep();
    }
}

void IOptimizer::fetchGradient() {
  assert("to be implemented");
}

bool IOptimizer::converged() {
  assert("to be implemented");
  return false;
}
