//
// Created by Michael Heuer on 26.01.17.
//


#include <iostream>
#include <iomanip>

#include "ElectronicWaveFunction.h"
#include "cppoptlib/solver/gradientdescentsolver.h"

class MyProblem : public cppoptlib::Problem<double,Eigen::Dynamic> {
public:
  //using typename cppoptlib::Problem<double, 18*3>::MyVector; // Inherit the Vector typedef

  MyProblem() {
    wf.setRandomElectronPositionCollection(18, ElectronPositioningMode::DENSITY);
    std::cout << "set random:" << std::endl;
    std::cout << wf.getElectronPositionCollection().transpose() << std::endl;
  }



  double value(const Eigen::VectorXd &x) {
    wf.evaluate(x);
    return wf.getNegativeLogarithmizedProbabilityDensity();
  }

  void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {

    wf.evaluate(x);
    grad = wf.getNegativeLogarithmizedSquaredElectronDriftCollection();
    //grad = wf.getElectronDriftCollection();
    assert(grad.cols() == 1);
    assert(x.cols() == 1);
    assert(grad.rows() == x.rows());

  }

  bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " ||dx|| = " << std::fixed << std::setw(8) << std::setprecision(4) << state.gradNorm
              << " ||x|| = "  << std::setw(6) << x.norm()
              << " f(x) = "   << std::setw(8) << value(x)
              << " x = [" << std::setprecision(8) << x.transpose() << "]" << std::endl;
    return true;
  }

  ElectronicWaveFunction wf;
};


int main(int argc, char const *argv[]) {

  Eigen::VectorXd  dummy(3*18);
  dummy << 0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, \
0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, \
0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, \
0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, \
0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53;

  MyProblem f;

  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults(); // Create a Criteria class to set the solver's stop conditions
  crit.iterations = 10000;                              // Increase the number of allowed iterations
  cppoptlib::GradientDescentSolver <MyProblem> solver;
  solver.setStopCriteria(crit);
  solver.minimize(f, dummy);
  std::cout << "f in argmin " << f(dummy) << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;

  return 0;
}

