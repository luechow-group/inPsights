//
// Created by Michael Heuer on 23.01.17.
//

#include <assert.h>
#include "Optimization/LBFGSOptimizer.h"

//#include "ElectronicWaveFunction.h"

LBFGSOptimizer::LBFGSOptimizer()
  : NewtonTypeOptimizer() {

  //wf_.createRandomElectronPositionCollection(18,ElectronPositioningMode::DENSITY); //TODO allow user to specify EPC

  // convert to positions vector
  //auto temp = wf_.getElectronPositionCollection();
  /*positions_ = Eigen::Map<Eigen::VectorXd>(temp.data(), temp.cols()*temp.rows());

  AX_.setcontent(positions_.size(), positions_.data());

  double epsg = 0.0000000001;
  double epsf = 0;
  double epsx = 0;
  alglib::ae_int_t maxits = 0;

  alglib::minlbfgscreate(positions_.size(), AX_, state_);
  alglib::minlbfgssetcond(state_, epsg, epsf, epsx, maxits);*/
  //alglib::minlbfgssetstpmax();
}

/*
void functionAndGradient(const alglib::real_1d_array &x,
                         double &func,
                         alglib::real_1d_array &grad, void *ptr)
{
}*/

void LBFGSOptimizer::performStep() {
  /*alglib::minlbfgsoptimize(state_, functionAndGradient);
  alglib::minlbfgsresults(state_, AX_, rep_);

  positions_ = Eigen::VectorXd(AX_.getcontent(),positions_.size());
   */
}


/*void LBFGSOptimizer::functionAndGradient(const alglib::real_1d_array &x,
                                         double &func,
                                         alglib::real_1d_array &grad, void *ptr)
{
  wf_.evaluate(positions_);
  func = wf_.getProbabilityDensity();
  auto temp = wf_.getElectronDriftCollection() * wf_.getProbabilityAmplitude();
  grad = temp.data();

}*/

void LBFGSOptimizer::constructHessian() {
  assert(false && "To be implemented");
}
