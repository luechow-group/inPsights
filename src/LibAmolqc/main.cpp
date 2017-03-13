#include <iostream>
#include <iomanip>
#include <ElectronicWaveFunction.h>
#include "problem.h"
#include "solver/bfgsnssolver.h"


class ElectronicWaveFunctionProblem : public cppoptlib::Problem<double,Eigen::Dynamic> {
public:

  ElectronicWaveFunctionProblem(){}

  double value(const Eigen::VectorXd &x) {
    wf.evaluate(x);
    return -wf.getProbabilityDensity();
    //return wf.getNegativeLogarithmizedProbabilityDensity();
  }

  void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    wf.evaluate(x);
    grad = -wf.getProbabilityDensityGradientCollection();
    //grad = wf.getNegativeLogarithmizedProbabilityDensityGradientCollection();
  }

  bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    Eigen::VectorXd grad;
    gradient(x,grad);
    std::cout << "(" << std::setw(2) << state.iterations << ")"
            << " f(x) = "     << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
            << " gradNorm = " << std::setw(8) << state.gradNorm
            //<< " xDelta = "   << std::setw(8) << state.xDelta
            << " g = [" << std::setprecision(16) << grad.transpose() << "]"
            //<< " x = [" << std::setprecision(16) << x.transpose() << "]"
            << std::endl;
  return true;
  }

private:
  ElectronicWaveFunction wf;
};

int main(int argc, char const *argv[]) {

  ElectronicWaveFunctionProblem f;


  Eigen::VectorXd x (2*3);
  x << 0.0, 0.0, 0.70, 0.0, 0.0, -0.70; //numerisch: klappt, analytisch: gradient läuft "schief"
  //x << 0.0, 0.0, 1.000, 0.0, 0.0, -0.270; //numerisch + analytisch: gradient aus convex hull search wird NAN

  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
  crit.iterations = 200;
  crit.gradNorm = 1e-6;
  cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
  solver.setDebug(cppoptlib::DebugLevel::High);
  solver.setStopCriteria(crit);
  solver.minimize(f, x);

  std::cout << "f in argmin " << f(x) << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;

  // f in argmin -0.0644447530112681 // with analytical gradient
  return 0;
}
