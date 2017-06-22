//
// Created by heuer on 06.04.17.
//

#include "ElectronicWaveFunctionProblem.h"

ElectronicWaveFunctionProblem::ElectronicWaveFunctionProblem()
        : wf_(ElectronicWaveFunction::getInstance()),valueCallCount_(0), gradientCallCount_(0) {}

double ElectronicWaveFunctionProblem::value(const Eigen::VectorXd &x) {
    valueCallCount_++;
    wf_.evaluate(x);
    return wf_.getNegativeLogarithmizedProbabilityDensity();
}

void ElectronicWaveFunctionProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    gradientCallCount_++;
    wf_.evaluate(x);
    grad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
}

void
ElectronicWaveFunctionProblem::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) {
  unsigned dims = x.size();

  std::vector<unsigned> ignoredIdx({0,1,2,3,4,9,10,11,12,13});
  assert(ignoredIdx.size() < dims);

  cppoptlib::Problem<double,Eigen::Dynamic>::hessian(x,hessian);

  std::cout << "reduced hessian" << std::endl;

  for (std::vector<unsigned>::const_iterator it = ignoredIdx.begin(); it != ignoredIdx.end(); ++it){
    assert( *it < dims);
    //set three rows to zero
    hessian.block(*it*3,0,3,dims) = Eigen::MatrixXd::Zero(3,dims);
    //set three columns to zero
    hessian.block(0,*it*3,dims,3) = Eigen::MatrixXd::Zero(dims,3);
  }
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
