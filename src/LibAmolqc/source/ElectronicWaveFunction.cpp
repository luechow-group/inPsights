//
// Created by Michael Heuer on 13.01.17.
//

#include "ElectronicWaveFunction.h"
#include <iostream>

extern "C" {
  void amolqc_init();
  void amolqc_set_wf(int *nElecs, int *nAtoms);
  void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
  void amolqc_set_electron_positions(double x[], int n);
  void amolqc_get_wave_function_values(double *phi, double *u);
  void amolqc_get_local_energy(double *elocal);
  void amolqc_get_drift(double drift[],int n);
}

ElectronicWaveFunction::ElectronicWaveFunction() {
  initialize();
}

ElectronicWaveFunction::~ElectronicWaveFunction() {
}

void ElectronicWaveFunction::setRandomElectronPositionCollection(unsigned electronNumber,
                                                                 ElectronPositioningMode::electronPositioningModeType
                                                                 electronPositioningModeType) {

  double *electronPositionCollectionArray = new double[electronNumber*3];
  amolqc_initial_positions(electronPositioningModeType, electronNumber,electronPositionCollectionArray);

  electronPositionCollection_ = Eigen::Map<Eigen::VectorXd>(electronPositionCollectionArray,
                                                            electronNumber*3, 1);

  delete electronPositionCollectionArray;
}

void ElectronicWaveFunction::initialize() {
  atomNumber_=0;
  electronNumber_=0;
  amolqc_init();
  amolqc_set_wf((int*)&electronNumber_, (int*)&atomNumber_);
  std::cout << electronNumber_ << ", " << atomNumber_ << std::endl;
  setRandomElectronPositionCollection(electronNumber_, ElectronPositioningMode::DENSITY);
}

double ElectronicWaveFunction::getLocalEnergy(){
  amolqc_get_local_energy(&localEnergy_);
  return localEnergy_;
}

void ElectronicWaveFunction::setElectronPositionCollection(const Eigen::VectorXd &electronPositionCollection) {
  assert(electronPositionCollection.rows() > 0);
  assert(electronPositionCollection.rows() % 3 == 0);

  electronPositionCollection_ = electronPositionCollection;
  Eigen::VectorXd copy = electronPositionCollection; // TODO ugly - redesign

  amolqc_set_electron_positions(copy.data(),electronNumber_);
};

void ElectronicWaveFunction::calculateWaveFunctionValues() {
  amolqc_get_wave_function_values(&determinantProbabilityAmplitude_, &jastrowFactor_);
}

void ElectronicWaveFunction::calculateWaveFunctionDrift() {

  double *electronDriftCollectionArray = new double[electronPositionCollection_.rows()];

  amolqc_get_drift(electronDriftCollectionArray, electronNumber_);

  electronDriftCollection_ = Eigen::Map<Eigen::VectorXd>( electronDriftCollectionArray,
                                                          electronPositionCollection_.rows(), 1);
  delete electronDriftCollectionArray;
}

double ElectronicWaveFunction::getDeterminantProbabilityAmplitude(){
  return determinantProbabilityAmplitude_;
};

double ElectronicWaveFunction::getJastrowFactor(){
  return jastrowFactor_;
};

double ElectronicWaveFunction::getProbabilityAmplitude(){
  return determinantProbabilityAmplitude_ * exp(jastrowFactor_);
};

double ElectronicWaveFunction::getProbabilityDensity(){
  return pow(getProbabilityAmplitude(),2);
};

double ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensity(){
  return -std::log(pow(getProbabilityAmplitude(),2));
};

Eigen::VectorXd ElectronicWaveFunction::getElectronPositionCollection(){
  return electronPositionCollection_;
};

Eigen::VectorXd ElectronicWaveFunction::getElectronDriftCollection(){
  return electronDriftCollection_;
};

Eigen::VectorXd ElectronicWaveFunction::getProbabilityAmplitudeGradientCollection(){
  return getProbabilityAmplitude()*electronDriftCollection_;
};

Eigen::VectorXd ElectronicWaveFunction::getProbabilityDensityGradientCollection(){
  return 2.0*getProbabilityAmplitude()* getProbabilityAmplitudeGradientCollection();
};
