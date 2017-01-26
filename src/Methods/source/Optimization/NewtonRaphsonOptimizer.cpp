//
// Created by Michael Heuer on 26.01.17.
//

#include "Optimization/NewtonRaphsonOptimizer.h"

NewtonRaphsonOptimizer::NewtonRaphsonOptimizer() {

}

void NewtonRaphsonOptimizer::performStep() {

  assert(false && "Hessian not calculated");
  assert(hessian_.rows() == hessian_.cols());
  assert(hessian_.rows() == gradients_.size());

  positions_ -= hessian_.inverse().eval() * gradients_;
}

void NewtonRaphsonOptimizer::constructHessian() {

  assert(false && "To be implemented");
}
