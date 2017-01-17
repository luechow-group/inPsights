//
// Created by Michael Heuer on 13.01.17.
//

#ifndef AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
#define AMOLQCGUI_ELECTRONICWAVEFUNCTION_H

#include <Eigen/Core>
//#include <cmath>

namespace ElectronPositioningMode {
  typedef enum {
    GAUSSIAN, DENSITY, LMO
  } electronPositioningModeType;
}

class ElectronicWaveFunction {
public:
  ElectronicWaveFunction();

  ~ElectronicWaveFunction();

  void initialize();

  void createRandomElectronPositionCollection(unsigned electronNumber,
                                              ElectronPositioningMode::electronPositioningModeType);

  void evaluate(Eigen::MatrixXd electronPositionCollection);

  double getLocalEnergy(){ return localEnergy_; };

  double getDeterminantProbabilityAmplitude(){ return determinantProbabilityAmplitude_; };

  double getJastrowFactor(){ return jastrowFactor_; };

  double getProbabilityAmplitude(){ return determinantProbabilityAmplitude_*jastrowFactor_; };

  double getProbabilityDensity(){ return pow(getProbabilityAmplitude(),2); };

  Eigen::MatrixXd getElectronDriftCollection(){ return electronDriftCollection_; };

  Eigen::MatrixXd getElectronPositionCollection(){ return electronPositionCollection_; };

private:
  void eigenMatrixToArray(Eigen::MatrixXd mat, double *arr);
  Eigen::MatrixXd arrayToEigenMatrix(double *arr, unsigned electronNumber);

  unsigned atomNumber_,electronNumber_;
  double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
  Eigen::MatrixXd electronPositionCollection_, electronDriftCollection_;

};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
