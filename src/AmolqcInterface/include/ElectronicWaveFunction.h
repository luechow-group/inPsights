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

extern "C" {
void amolqc_init();
void amolqc_set_wf(int *nElecs, int *nAtoms);
void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
void amolqc_eloc(double x[], int n, double *phi, double *u, double grad[], double *elocal);
}

class ElectronicWaveFunction {

public:
  static ElectronicWaveFunction& getInstance();

  void initialize();

  void setRandomElectronPositionCollection(unsigned electronNumber,
                                           ElectronPositioningMode::electronPositioningModeType);

  void evaluate(const Eigen::VectorXd &electronPositionCollection);

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

  Eigen::VectorXd getNegativeLogarithmizedProbabilityDensityGradientCollection();

private:
  ElectronicWaveFunction();
  unsigned atomNumber_,electronNumber_;
  double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
  Eigen::VectorXd electronPositionCollection_, electronDriftCollection_;

};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
