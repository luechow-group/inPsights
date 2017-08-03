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
void amolqc_set_wf(int *nElecs, int *nAtoms, const char* fileName);
void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
void amolqc_eloc(double x[], int n, double *phi, double *u, double grad[], double *elocal);
}

class ElectronicWaveFunction {

public:
  static ElectronicWaveFunction& getInstance(const std::string& fileName = "");

  const std::string& getFileName();

  void initialize(const std::string& fileName);

  void setRandomElectronPositionCollection(unsigned electronNumber,
                                           ElectronPositioningMode::electronPositioningModeType);

  void evaluate(const Eigen::VectorXd &electronPositionCollection);

  double getLocalEnergy();

  double getDeterminantProbabilityAmplitude();

  double getJastrowFactor();

  double getProbabilityAmplitude();

  double getProbabilityDensity();

  double getNegativeLogarithmizedProbabilityDensity();

    double getInverseNegativeLogarithmizedProbabilityDensity();

  Eigen::VectorXd getElectronPositionCollection();

  Eigen::VectorXd getElectronDriftCollection();

  Eigen::VectorXd getProbabilityAmplitudeGradientCollection();

  Eigen::VectorXd getProbabilityDensityGradientCollection();

  Eigen::VectorXd getNegativeLogarithmizedProbabilityDensityGradientCollection();

    Eigen::VectorXd getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();

  int getNumberOfNuclei() const;

  int getNumberOfElectrons() const;

  void setFrozenElectrons(const std::vector<int>& frozenElectrons);

private:
  ElectronicWaveFunction(const std::string& fileName);
  const std::string fileName_;
  unsigned numberOfNuclei_, numberOfElectrons_;
  double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
  Eigen::VectorXd electronPositionCollection_, electronDriftCollection_;

  std::vector<int> frozenElectrons_;
};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
