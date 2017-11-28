#include <iostream>
#include <iomanip>
#include "ElectronicWaveFunctionProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/timeintegrationsolver.h"
#include <Eigen/Eigenvalues>

int main(int argc, char const *argv[]) {

  ElectronicWaveFunction::getInstance("Ethane-em-5.wf");

  /*amolqcInput1*/
  /*Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xA << \
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
  0.005195,  0.207915, -1.906905;*/

  /*amolqcOutput1 << \
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
  0.010891, -0.034989, -0.645852;*/

  //Ethane global max
  Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xA <<
     0.000000, 0.000000, 1.443184,\
     0.000000, 0.000000,-1.443184,\
    -1.662146,-0.959641, 2.192989,\
    -1.662146, 0.959641,-2.192989,\
     0.000000,-1.919300,-2.192989,\
    -0.024099, 0.773535, 1.718336+0.3,\
     0.657845,-0.407636, 1.718335-0.3,\
     0.658424, 0.380140,-1.752435+0.3,\
    -0.034848,-0.020119,-0.660100-0.3,\
     0.000000, 0.000000, 1.443184,\
     0.000000, 0.000000,-1.443184,\
     1.662146,-0.959641, 2.192989,\
     1.662146, 0.959641,-2.192989,\
     0.000000, 1.919300, 2.192989,\
    -0.657846, 0.407635,-1.718336-0.3,\
     0.024100,-0.773537,-1.718336+0.3,\
    -0.658423,-0.380140, 1.752435-0.3,\
     0.034847, 0.020120, 0.660096+0.3;


  ElectronicWaveFunctionProblem f("");
  Eigen::VectorXd grad(18*3);
  grad.setZero(18*3);
  f.gradient(xA,grad);
  std::cout << std::setprecision(1) << grad.transpose() << std::endl;

  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
  //crit.iterations = 1000;
  crit.gradNorm = 1e-10;
  //cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
  cppoptlib::TimeIntegrationSolver<ElectronicWaveFunctionProblem> solver;
  solver.setDebug(cppoptlib::DebugLevel::High);
  solver.setStopCriteria(crit);

  Eigen::VectorXd x = xA;
  solver.minimize(f, x);


  /*std::cout << "f in argmin " << f(x) << std::endl;
  std::cout << "Solver status: " << solver.status() << std::endl;
  std::cout << "Final criteria values: " << std::endl << solver.criteria() << std::endl;*/
  std::cout << "Total eloc calls: " << f.getTotalElocCalls() << std::endl;

  grad.setZero(18*3);
  f.gradient(xA,grad);
  std::cout << std::setprecision(1) << grad.transpose() << std::endl;

  Eigen::MatrixXd hess(18*3,18*3);
  hess.setZero(18*3,18*3);
  f.hessian(xA,hess);
  std::cout << std::setprecision(1) << hess << std::endl;

  Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess,false);
  auto eigenvalues = eigenSolver.eigenvalues();
  std::cout << eigenvalues << std::endl;

  return 0;
}
