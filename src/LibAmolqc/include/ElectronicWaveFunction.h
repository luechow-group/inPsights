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
  double getLocalEnergy(){ return localEnergy; };
  double getDeterminantProbabilityAmplitude(){ return determinantProbabilityAmplitude; };
  double getJastrowFactor(){ return jastrowFactor; };
  double getProbabilityAmplitude(){ return determinantProbabilityAmplitude*jastrowFactor; };
  double getProbabilityDensity(){ return pow(getProbabilityAmplitude(),2); };
  Eigen::MatrixXd getElectronDriftCollection(){ return electronDriftCollection; };
  Eigen::MatrixXd getElectronPositionCollection(){ return electronPositionCollection; };

private:
  unsigned atomNumber,electronNumber;
  double determinantProbabilityAmplitude, jastrowFactor, localEnergy;
  Eigen::MatrixXd electronPositionCollection, electronDriftCollection;//electronGradientCollection;

  void eigenMatrixToArray(Eigen::MatrixXd mat, double *arr);
  Eigen::MatrixXd arrayToEigenMatrix(double *arr, unsigned electronNumber);
};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
