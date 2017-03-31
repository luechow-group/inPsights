//
// Created by Michael Heuer on 13.01.17.
//

#include "ElectronicWaveFunction.h"
#include <iostream>

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

void ElectronicWaveFunction::evaluate(const Eigen::VectorXd &electronPositionCollection) {
  assert(electronPositionCollection.rows() > 0);
  assert(electronPositionCollection.rows() % 3 == 0);

  electronPositionCollection_ = electronPositionCollection;

  double *electronDriftCollectionArray = new double[electronPositionCollection.rows()];

  Eigen::VectorXd copy = electronPositionCollection; // TODO ugly - redesign

  amolqc_eloc(copy.data(), electronNumber_, &determinantProbabilityAmplitude_, &jastrowFactor_,
              electronDriftCollectionArray, &localEnergy_);

  electronDriftCollection_ = Eigen::Map<Eigen::VectorXd>( electronDriftCollectionArray,
                                                          electronPositionCollection.rows(), 1);
  delete electronDriftCollectionArray;
}

double ElectronicWaveFunction::getLocalEnergy(){
  return localEnergy_;
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
  std::cout << electronDriftCollection_ << std::endl;
  return electronDriftCollection_;
};

Eigen::VectorXd ElectronicWaveFunction::getProbabilityAmplitudeGradientCollection(){
  return getProbabilityAmplitude()*electronDriftCollection_;
};

Eigen::VectorXd ElectronicWaveFunction::getProbabilityDensityGradientCollection(){
  return 2.0*getProbabilityAmplitude()* getProbabilityAmplitudeGradientCollection();
};

Eigen::VectorXd ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensityGradientCollection() {
  return -2.0 * electronDriftCollection_;
}
