//
// Created by Michael Heuer on 13.01.17.
//

#include "ElectronicWaveFunction.h"


/*TODO Remember to set the amolqc path for runtime in the configuration file
 ** AMOLQC=/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc/
 *TODO Remember to put the t.wf file in the executable folders
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-debug/src/LibAmolqc/t.wf
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-release/src/LibAmolqc/t.wf
 */
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

void ElectronicWaveFunction::createRandomElectronPositionCollection(unsigned electronNumber,
                                                                    ElectronPositioningMode::electronPositioningModeType
                                                                    electronPositioningModeType) {

  double *electronPositionCollectionArray = new double[electronNumber*3];
  amolqc_initial_positions(electronPositioningModeType, electronNumber,electronPositionCollectionArray);

  electronPositionCollection_ = arrayToEigenMatrix(electronPositionCollectionArray, electronNumber);

  delete electronPositionCollectionArray;
}

void ElectronicWaveFunction::initialize() {
  atomNumber_=0;
  electronNumber_=0;
  amolqc_init();
  amolqc_set_wf((int*)&electronNumber_, (int*)&atomNumber_);

  createRandomElectronPositionCollection(electronNumber_, ElectronPositioningMode::DENSITY);
}

void ElectronicWaveFunction::evaluate(Eigen::MatrixXd electronPositionCollection) {
  assert(electronPositionCollection.rows() > 0);
  assert(electronPositionCollection.cols() == 3);

  electronPositionCollection_ = electronPositionCollection;

  long electronCoordinatesNumber = electronPositionCollection.rows()* electronPositionCollection.cols();
  double *electronPositionCollectionArray = new double[electronCoordinatesNumber];
  double *electronDriftCollectionArray= new double[electronCoordinatesNumber];

  eigenMatrixToArray(electronPositionCollection, electronPositionCollectionArray);
  amolqc_eloc(electronPositionCollectionArray, electronNumber_, &determinantProbabilityAmplitude_, &jastrowFactor_, electronDriftCollectionArray, &localEnergy_);


  electronDriftCollection_ = arrayToEigenMatrix(electronDriftCollectionArray, electronNumber_);

  delete electronPositionCollectionArray;
  delete electronDriftCollectionArray;
}


void ElectronicWaveFunction::eigenMatrixToArray(Eigen::MatrixXd mat, double *arr) {
  for (unsigned i = 0; i < mat.rows(); ++i) {
    for (unsigned j = 0; j < mat.cols(); ++j) {
      arr[i*3 + j] = mat(i,j);
    }
  }
  /*
  delete arr;
  arr = nullptr;
  Eigen::Map<Eigen::MatrixXd>( arr, mat.rows(), mat.cols() ) = mat;
   */
}

Eigen::MatrixXd ElectronicWaveFunction::arrayToEigenMatrix(double *arr, unsigned electronNumber) {
  /*Eigen::MatrixXd mat(electronNumber,3);
  for (unsigned i = 0; i < electronNumber; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      mat(i, j) = arr[i * 3 + j];
    }
  }*/

  return Eigen::Map<Eigen::MatrixXd>(arr, electronNumber, 3);
}


