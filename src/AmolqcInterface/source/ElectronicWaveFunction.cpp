//
// Created by Michael Heuer on 13.01.17.
//

#include "ElectronicWaveFunction.h"
#include <iostream>

ElectronicWaveFunction &ElectronicWaveFunction::getInstance(const std::string& fileName) {

  // these members are static and thus only initialized once
  static ElectronicWaveFunction electronicWaveFunction(fileName);

  if( fileName != electronicWaveFunction.getFileName() && fileName != "" )
    std::cout << "The current wavefunction file is " << electronicWaveFunction.getFileName() << ".\n"
              << " It cannot be reinitialized to " << fileName << "."
              << std::endl;

  return electronicWaveFunction;
}

const std::string &ElectronicWaveFunction::getFileName() {
  return fileName_;
}

ElectronicWaveFunction::ElectronicWaveFunction(const std::string& fileName)
        :fileName_(fileName) {
  initialize(fileName);
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

void ElectronicWaveFunction::initialize(const std::string& fileName) {
  numberOfNuclei_=0;
  numberOfElectrons_=0;

  char* f = (char*) fileName.c_str();
  std::cout << f << std::endl;

  amolqc_init();
  amolqc_set_wf((int*)&numberOfElectrons_, (int*)&numberOfNuclei_, f);
  std::cout << numberOfElectrons_ << ", " << numberOfNuclei_ << std::endl;
  setRandomElectronPositionCollection(numberOfElectrons_, ElectronPositioningMode::DENSITY);
}

int ElectronicWaveFunction::getNumberOfNuclei() const {
  return numberOfNuclei_;
}

int ElectronicWaveFunction::getNumberOfElectrons() const {
  return numberOfElectrons_;
}

void ElectronicWaveFunction::evaluate(const Eigen::VectorXd &electronPositionCollection) {
  assert(electronPositionCollection.rows() > 0);
  assert(electronPositionCollection.rows() % 3 == 0);

  electronPositionCollection_ = electronPositionCollection;

  double *electronDriftCollectionArray = new double[electronPositionCollection.rows()];

  Eigen::VectorXd copy = electronPositionCollection; // TODO ugly - redesign

  amolqc_eloc(copy.data(), numberOfElectrons_, &determinantProbabilityAmplitude_, &jastrowFactor_,
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
  return determinantProbabilityAmplitude_ * std::exp(jastrowFactor_); ///
  //std::exp(-339); //ethane
  //std::exp(-11); // ethylene
};

double ElectronicWaveFunction::getProbabilityDensity(){
  return pow(getProbabilityAmplitude(),2);
};

double ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensity(){
  return -std::log(pow(getProbabilityAmplitude(),2));
};

double ElectronicWaveFunction::getInverseNegativeLogarithmizedProbabilityDensity() {
  //return 1.0/ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensity();
  return -1.0/std::log(pow(getProbabilityAmplitude(),2));
}

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
  for (std::vector<int>::const_iterator it = frozenElectrons_.begin(); it != frozenElectrons_.end(); ++it ){
    // index of the elctrons start at 0
    electronDriftCollection_.segment(((*it)-1)*3,3) = Eigen::VectorXd::Zero(3);
  }

  return 2.0*getProbabilityAmplitude()* getProbabilityAmplitudeGradientCollection();
};

Eigen::VectorXd ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensityGradientCollection() {
  for (std::vector<int>::const_iterator it = frozenElectrons_.begin(); it != frozenElectrons_.end(); ++it ){
    // index of the elctrons start at 0
    electronDriftCollection_.segment(((*it)-1)*3,3) = Eigen::VectorXd::Zero(3);
  }
  return -2.0 * electronDriftCollection_;

}

Eigen::VectorXd ElectronicWaveFunction::getInverseNegativeLogarithmizedProbabilityDensityGradientCollection() {
  //return getNegativeLogarithmizedProbabilityDensityGradientCollection();//.cwiseInverse();
  return (-2.0 * electronDriftCollection_).cwiseInverse();
  //return 0.5*getProbabilityAmplitude()*getProbabilityAmplitudeGradientCollection().cwiseInverse();
}

void ElectronicWaveFunction::setFrozenElectrons(const std::vector<int>& frozenElectrons) {
  frozenElectrons_= frozenElectrons;
}

std::vector<int> ElectronicWaveFunction::getFrozenElectrons() {
  return frozenElectrons_;
}
