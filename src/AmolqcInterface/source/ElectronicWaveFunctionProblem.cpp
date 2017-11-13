//
// Created by heuer on 06.04.17.
//

#include "ElectronicWaveFunctionProblem.h"

ElectronicWaveFunctionProblem::ElectronicWaveFunctionProblem(const std::string fileName)
        : wf_(ElectronicWaveFunction::getInstance(fileName)),valueCallCount_(0), gradientCallCount_(0) {}

double ElectronicWaveFunctionProblem::value(const Eigen::VectorXd &x) {
    valueCallCount_++;
    wf_.evaluate(x);
  return wf_.getNegativeLogarithmizedProbabilityDensity();
  //return wf_.getProbabilityDensity();

}

void ElectronicWaveFunctionProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  gradientCallCount_++;
  wf_.evaluate(x);
  grad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
  //grad = wf_.getProbabilityDensityGradientCollection();
}

void ElectronicWaveFunctionProblem::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) {
  unsigned dims = x.size();


  std::vector<int> ignoredIdx = wf_.getFrozenElectrons();
  assert(ignoredIdx.size() < dims);

  cppoptlib::Problem<double,Eigen::Dynamic>::hessian(x,hessian);

  //std::cout << "orig hessian " << hessian << std::endl;

  for (std::vector<int>::const_iterator it = ignoredIdx.begin(); it != ignoredIdx.end(); ++it){
    assert( *it < dims);
    //set three rows to zero
    hessian.block(*it*3,0,3,dims) = Eigen::MatrixXd::Zero(3,dims);
    //set three columns to zero
    hessian.block(0,*it*3,dims,3) = Eigen::MatrixXd::Zero(dims,3);
  }
  //std::cout << "reduced hessian " << hessian << std::endl;
}

bool ElectronicWaveFunctionProblem::callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    //notifyObserversAboutPerformedStep();

    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " f(x) = "     << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
              << " xDelta = "   << std::setw(8) << state.xDelta
              << " gradInfNorm = "   << std::setw(8) << state.gradNorm
              << std::endl;
    return true;
}
