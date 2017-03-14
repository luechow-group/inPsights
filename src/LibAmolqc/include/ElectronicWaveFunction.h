//
// Created by Michael Heuer on 13.01.17.
//

#ifndef AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
#define AMOLQCGUI_ELECTRONICWAVEFUNCTION_H

#include <Eigen/Core>

namespace ElectronPositioningMode {
  typedef enum {
    GAUSSIAN, DENSITY, LMO
  } electronPositioningModeType;
}

/*TODO Remember to set the amolqc path for runtime in the configuration file
 ** AMOLQC=/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc/
 *TODO Remember to put the t.wf file in the executable folders
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-debug/src/LibAmolqc/t.wf
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-release/src/LibAmolqc/t.wf
 */


class ElectronicWaveFunction {

public:
  ElectronicWaveFunction();

  ~ElectronicWaveFunction();

  void initialize();

  void setRandomElectronPositionCollection(unsigned electronNumber,
                                           ElectronPositioningMode::electronPositioningModeType);

  void setElectronPositionCollection(const Eigen::VectorXd &electronPositionCollection);

  void calculateWaveFunctionValues();

  void calculateWaveFunctionDrift();

  double getLocalEnergy();

  double getDeterminantProbabilityAmplitude();

  double getJastrowFactor();

  double getProbabilityAmplitude();

  double getProbabilityDensity();

  double getNegativeLogarithmizedProbabilityDensity();

  Eigen::VectorXd getElectronPositionCollection();

  Eigen::VectorXd getElectronDriftCollection();

  Eigen::VectorXd getProbabilityAmplitudeGradientCollection();

  Eigen::VectorXd getProbabilityDensityGradientCollection();

  //Eigen::VectorXd getNegativeLogarithmizedProbabilityDensityGradientCollection();
  // not allowed - gradient can be negative

private:
  unsigned atomNumber_,electronNumber_;
  double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
  Eigen::VectorXd electronPositionCollection_, electronDriftCollection_;

};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
