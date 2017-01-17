#include <iostream>
#include <ElectronicWaveFunction.h>

int main(int argc, char const *argv[]) {

  ElectronicWaveFunction wf;

  Eigen::MatrixXd test (18,3);
  //test.setRandom();
  test << 0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, \
0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, \
0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, \
0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, \
0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53;

  //wf.createRandomElectronPositionCollection(18,ElectronPositioningMode::DENSITY);
  //wf.evaluate(wf.getElectronPositionCollection());
  wf.evaluate(test);
  std::cout << wf.getElectronPositionCollection() << std::endl << std::endl;
  std::cout << wf.getLocalEnergy() << std::endl << std::endl;
  std::cout << wf.getDeterminantProbabilityAmplitude() << std::endl << std::endl;
  std::cout << wf.getJastrowFactor() << std::endl << std::endl;
  std::cout << wf.getProbabilityDensity() << std::endl << std::endl;
  std::cout << wf.getElectronDriftCollection() << std::endl << std::endl;
  return 0;
}

