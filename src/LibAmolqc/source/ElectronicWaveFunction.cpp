//
// Created by Michael Heuer on 13.01.17.
//

#include "ElectronicWaveFunction.h"
#include <iostream>

extern "C" {
  void amolqc_init();
  void amolqc_set_wf(int *nElecs, int *nAtoms);
  void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
  void amolqc_eloc(double x[], int n, double *phi, double *u, double grad[], double *elocal);
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

/*void ElectronicWaveFunction::initialize(const Eigen::VectorXd &electronPositionCollection) {
  atomNumber_=0;
  electronNumber_=0;
  amolqc_init();
  amolqc_set_wf((int*)&electronNumber_, (int*)&atomNumber_);

  setElectronPositionCollection(electronPositionCollection); //TODO IMPLEMENT
}*/

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


