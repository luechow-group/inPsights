#include <iostream>
#include <iomanip>
#include <ElectronicWaveFunction.h>
#include "problem.h"
#include "solver/bfgsnssolver.h"


class ElectronicWaveFunctionProblem : public cppoptlib::Problem<double,Eigen::Dynamic> {
public:

  ElectronicWaveFunctionProblem()
    : valueCallCount_(0), gradientCallCount_(0) {}

  double value(const Eigen::VectorXd &x) {
    valueCallCount_++;
    wf_.evaluate(x);
    return wf_.getNegativeLogarithmizedProbabilityDensity();
  }

  void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    gradientCallCount_++;
    wf_.evaluate(x);
    grad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
  }

  bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    Eigen::VectorXd grad;
    gradient(x,grad);
    std::cout << "(" << std::setw(2) << state.iterations << ")"
            << " f(x) = "     << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
            << " gradNorm = " << std::setw(8) << state.gradNorm
            << " xDelta = "   << std::setw(8) << state.xDelta
            //<< " g = [" << std::setprecision(16) << grad.transpose() << "]"
            //<< " x = [" << std::setprecision(16) << x.transpose() << "]"
            << std::endl;
  return true;
  }

  unsigned getValueCallCount(){
    return valueCallCount_;
  }

  unsigned getGradientCallCount(){
    return gradientCallCount_;
  }

  unsigned getTotalElocCalls(){
    return getValueCallCount()+getGradientCallCount();
  }

  void resetCounters(){
    valueCallCount_ = 0;
    gradientCallCount_ = 0;
  }

private:
  unsigned valueCallCount_, gradientCallCount_;
  ElectronicWaveFunction wf_;
};

int main(int argc, char const *argv[]) {

  Eigen::VectorXd amolqcInput1(18*3);
  Eigen::VectorXd amolqcOutput1(18*3);

  amolqcInput1 << \
  0.714583,  2.171709,  2.377429,\
 -0.805267,  0.373607,  0.961730,\
 -0.013201, -0.133104,  1.591434,\
  0.305421,  0.141948, -1.348573,\
  0.733193,  0.981820, -1.175646,\
 -1.631408,  0.621465, -2.145863,\
 -0.549092, -2.641827,  3.075085,\
  1.668276,  1.311450, -0.564745,\
  1.554871, -1.617958, -2.978076,\
  0.204682, -0.074675,  1.290531,\
 -0.880830, -0.567550,  0.091367,\
  0.483278, -2.256104,  1.174406,\
  1.888764, -0.491579, -1.046459,\
 -0.154281,  1.014234, -2.217571,\
 -0.924312,  0.945934, -0.019794,\
 -1.987497,  0.072370,  1.736939,\
 -0.544636, -2.204059, -3.499582,\
  0.005195,  0.207915, -1.906905;

  amolqcOutput1 << \
  0.000000,  0.000000,  1.446226,\
  0.000000,  0.000000, -1.446226,\
  0.961606,  1.667381,  2.199093,\
 -0.961606,  1.667381, -2.199093,\
  0.963193, -1.666474,  2.199093,\
 -0.435950, -0.706833, -1.719891,\
 -0.772256, -0.007354,  1.744558,\
  0.034146, -0.001214,  0.706588,\
  0.829139,  0.022591, -1.725390,\
  0.000000,  0.000000,  1.446226,\
  0.000000,  0.000000, -1.446226,\
  0.963193, -1.666474,  2.199093,\
 -0.963193, -1.666474, -2.199093,\
  1.924799, -0.000907, -2.199093,\
 -1.924799, -0.000907,  2.199093,\
  0.394929,  0.672687,  1.724599,\
 -0.414574,  0.718595, -1.753196,\
  0.010891, -0.034989, -0.645852;


  ElectronicWaveFunctionProblem f;

  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
  crit.iterations = 1000;
  crit.xDelta = 1e-4; // for 1e-5 the gradient in convex hull opt becomes NaN
  cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
  solver.setDebug(cppoptlib::DebugLevel::High);
  solver.setStopCriteria(crit);

  Eigen::VectorXd x = amolqcInput1;
  solver.minimize(f, x);

  // Permute results
  Eigen::VectorXi permutation1(18);
  permutation1 << 3, 4, 1, 5, 7, 6, 2, 8, 9, 10, 18, 12, 17, 13, 16, 15, 14, 11;
  // correct for zero-based counting
  permutation1.array() -= 1;

  Eigen::VectorXd optimResultPermuted1(18*3);

  for (int i = 0; i < 18 ; ++i) {
    optimResultPermuted1.segment(i*3,3) = x.segment(permutation1(i)*3,3);
  }

  std::cout << std::setprecision(7);
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>> (amolqcOutput1.data()) << std::endl;
  std::cout << std::endl;
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>> (optimResultPermuted1.data()) << std::endl;
  std::cout << std::endl;
  std::cout << "diff to amolqc result:" << std::endl;
  Eigen::VectorXd diff(optimResultPermuted1 - amolqcOutput1);
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>>(diff.data())<< std::endl;
  std::cout << std::endl;
  std::cout << "f in argmin " << f(x) << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;
  std::cout << "Total eloc calls: " << f.getTotalElocCalls() << std::endl;

  /*
   *
   */

  Eigen::VectorXd amolqcInput2(18*3);
  Eigen::VectorXd amolqcOutput2(18*3);

  amolqcInput2 << \
  0.925661,  1.099722,  1.941702,\
 -0.550422,  0.237153,  1.118380,\
  0.023390,  0.218214,  1.480570,\
 -0.092484,  0.019468, -1.352853,\
  0.813849,  0.298597, -1.738673,\
 -0.706151, -0.087879, -2.678484,\
 -1.192049, -1.511869,  2.971718,\
  1.838029,  0.898649,  0.598746,\
  1.694989, -1.735890, -3.939891,\
  0.098824, -0.072349,  1.701301,\
 -0.899394, -0.369229, -0.123691,\
  1.136408, -1.748549,  0.877996,\
  0.613289, -0.220153, -2.066445,\
 -0.785748,  1.522570, -2.252711,\
 -0.385672,  1.823523,  0.615760,\
 -2.044852, -0.187601,  1.530266,\
 -0.981077, -2.158334, -3.623649,\
 -0.100396,  0.127356, -1.320765;

  amolqcOutput2 << \
  0.000000,  0.000000,  1.446226,\
  0.000000,  0.000000, -1.446226,\
  0.961606,  1.667381,  2.199093,\
  1.924799, -0.000907, -2.199093,\
 -0.406422,  0.746455, -1.739065,\
 -0.407068, -0.747451, -1.733527,\
  0.439809, -0.723812,  1.760897,\
  0.025154,  0.028000,  0.672916,\
 -0.846128,  0.018941,  1.754514,\
  0.000000,  0.000000,  1.446226,\
  0.000000,  0.000000, -1.446226,\
 -1.924799, -0.000907,  2.199093,\
  0.963193, -1.666474,  2.199093,\
 -0.961606,  1.667381, -2.199093,\
 -0.963193, -1.666474, -2.199093,\
  0.840611, -0.000591, -1.764839,\
  0.423257,  0.734397,  1.746271,\
  0.039757, -0.008220, -0.667437;


  cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver2;
  solver2.setDebug(cppoptlib::DebugLevel::High);
  solver2.setStopCriteria(crit);

  x = amolqcInput2;
  f.resetCounters();
  solver2.minimize(f, x);

  Eigen::VectorXi permutation2(18);
  permutation2 << 3, 4, 1, 9, 5, 6, 7, 8, 2, 10, 18, 16, 12, 14, 17, 13, 15, 11;
  // correct for zero-based counting
  permutation2.array() -= 1;

  Eigen::VectorXd optimResultPermuted2(18*3);

  for (int i = 0; i < 18 ; ++i) {
    optimResultPermuted2.segment(i*3,3) = x.segment(permutation2(i)*3,3);
  }

  std::cout << std::setprecision(7);
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>> (amolqcOutput2.data()) << std::endl;
  std::cout << std::endl;
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>> (optimResultPermuted2.data()) << std::endl;
  std::cout << std::endl;
  std::cout << "diff to amolqc result:" << std::endl;
  Eigen::VectorXd diff2(optimResultPermuted2 - amolqcOutput2);
  std::cout << Eigen::Map<Eigen::Matrix<double ,18,3,Eigen::RowMajor>>(diff2.data())<< std::endl;
  std::cout << std::endl;
  std::cout << "f in argmin " << f(x) << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;
  std::cout << "Total eloc calls: " << f.getTotalElocCalls() << std::endl;

  /*
  ElectronicWaveFunctionProblem f;
  Eigen::VectorXd x (2*3);
  x << 0.0, 1.0, 0.70, -2.0, 0.0, -0.70;

  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
  crit.iterations = 200;
  //crit.xDelta = 1e-6;
  cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
  solver.setDebug(cppoptlib::DebugLevel::High);
  solver.setStopCriteria(crit);
  solver.minimize(f, x);

  std::cout << "f in argmin " << f(x) << std::endl;
  std::cout << " x = [" << std::setprecision(16) << x.transpose() << "]" << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;
*/

  return 0;
}
