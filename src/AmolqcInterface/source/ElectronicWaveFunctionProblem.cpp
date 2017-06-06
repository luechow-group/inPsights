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

bool ElectronicWaveFunctionProblem::callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    notifyObserversAboutPerformedStep();

    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " f(x) = "     << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
              << " xDelta = "   << std::setw(8) << state.xDelta
              << std::endl;
    return true;
}
