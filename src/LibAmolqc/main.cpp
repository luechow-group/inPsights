#include <iostream>
#include <ElectronicWaveFunction.h>

int main(int argc, char const *argv[]) {

  ElectronicWaveFunction wf;

  Eigen::MatrixXd test (18,3);
  test.setRandom();
  //wf.createRandomElectronPositionCollection(18,ElectronPositioningMode::DENSITY);
  //wf.evaluate(wf.getElectronPositionCollection());
  wf.evaluate(test);
  std::cout << wf.getLocalEnergy() << std::endl;
  std::cout << wf.getDeterminantProbabilityAmplitude() << std::endl;
  std::cout << wf.getJastrowFactor() << std::endl;
  std::cout << wf.getProbabilityDensity() << std::endl;
  std::cout << wf.getElectronDriftCollection() << std::endl;
  return 0;
}

