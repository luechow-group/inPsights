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
  //void initialize(const Eigen::VectorXd &electronPositionCollection);

  void setRandomElectronPositionCollection(unsigned electronNumber,
                                           ElectronPositioningMode::electronPositioningModeType);

  void evaluate(const Eigen::VectorXd &electronPositionCollection);

  double getLocalEnergy(){ return localEnergy_; };

  double getDeterminantProbabilityAmplitude(){ return determinantProbabilityAmplitude_; };

  double getJastrowFactor(){ return jastrowFactor_; };

  double getProbabilityAmplitude(){ return determinantProbabilityAmplitude_ * exp(jastrowFactor_); };

  double getProbabilityDensity(){ return pow(getProbabilityAmplitude(),2); };

  double getNegativeLogarithmizedProbabilityDensity(){ return -log(pow(getProbabilityAmplitude(),2)); };

  Eigen::VectorXd getElectronPositionCollection(){ return electronPositionCollection_; };

  Eigen::VectorXd getElectronDriftCollection(){ return electronDriftCollection_; };

  Eigen::VectorXd getSquaredElectronDriftCollection(){ return electronDriftCollection_.array().square(); };

  Eigen::VectorXd getNegativeLogarithmizedSquaredElectronDriftCollection(){
    return -Eigen::log(getSquaredElectronDriftCollection().array());
  };



private:
  unsigned atomNumber_,electronNumber_;
  double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
  Eigen::VectorXd electronPositionCollection_, electronDriftCollection_;

};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
