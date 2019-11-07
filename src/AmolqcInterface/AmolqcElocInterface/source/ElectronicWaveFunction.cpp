/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include "ElectronicWaveFunction.h"
#include "WfFileImporter.h"

ElectronicWaveFunction &ElectronicWaveFunction::getInstance(const std::string& fileName) {

  // these members are static and thus only initialized once
  static ElectronicWaveFunction electronicWaveFunction(fileName);

  if( fileName != electronicWaveFunction.getFileName() && fileName != "" )
    std::cout << "The current wavefunction file is " << electronicWaveFunction.getFileName() << ".\n"
              << " It cannot be reinitialized to " << fileName << "."
              << std::endl;

  return electronicWaveFunction;
}

ElectronicWaveFunction &ElectronicWaveFunction::getEmpty() {

    static ElectronicWaveFunction electronicWaveFunction;

    return electronicWaveFunction;
}

const std::string &ElectronicWaveFunction::getFileName() {
  return fileName_;
}

ElectronicWaveFunction::ElectronicWaveFunction()
        :fileName_(""),
         numberOfNuclei_(0),
         numberOfElectrons_(0),
         numberOfAlphaElectrons_(0),
         numberOfBetaElectrons_(0),
         determinantProbabilityAmplitude_(0),
         jastrowFactor_(0),
         localEnergy_(0),
         electronPositionsVectorAsEigenVector_(0),
         electronDriftCollection_(0),
         atomsVector_(AtomsVector()),
         spinTypesVector_(0)
{
    amolqc_init();
}

ElectronicWaveFunction::ElectronicWaveFunction(const std::string& fileName)
        :fileName_(fileName) {
  initialize(fileName);
}

void ElectronicWaveFunction::setRandomElectronPositionsVector(unsigned electronNumber,
                                                                 ElectronPositioningMode::electronPositioningModeType
                                                                 electronPositioningModeType) {
  auto *electronPositionsVectorArray = new double[electronNumber*3];
  amolqc_initial_positions(electronPositioningModeType, electronNumber, electronPositionsVectorArray);

  electronPositionsVectorAsEigenVector_ = Eigen::Map<Eigen::VectorXd>(electronPositionsVectorArray,
                                                            electronNumber*3, 1);
  delete electronPositionsVectorArray;
}

void ElectronicWaveFunction::initialize(const std::string& fileName) {
  numberOfNuclei_=0;
  numberOfElectrons_=0;

  auto f = (char*) fileName.c_str();
  std::cout << f << std::endl;

  amolqc_init();
  amolqc_set_wf((int*)&numberOfElectrons_, (int*)&numberOfNuclei_, f);
  setRandomElectronPositionsVector((unsigned)numberOfElectrons_, ElectronPositioningMode::DENSITY);

  // import atoms from .wf file
  WfFileImporter wfFileImporter(fileName);
  assert(numberOfElectrons_ == wfFileImporter.getNumberOfElectrons()
         && "Number of electrons from amolqc must match the number by counting the nuclear charges minus the overall charge of the molecule");
  numberOfBetaElectrons_ = wfFileImporter.getNumberOfBetaElectrons();
  numberOfAlphaElectrons_ = wfFileImporter.getNumberOfAlphaElectrons();
  atomsVector_ = wfFileImporter.getAtomsVector();
  spinTypesVector_  = SpinTypesVector(numberOfAlphaElectrons_,numberOfBetaElectrons_);
}

unsigned long ElectronicWaveFunction::getNumberOfNuclei() const {
  return numberOfNuclei_;
}

unsigned long ElectronicWaveFunction::getNumberOfElectrons() const {
  return numberOfElectrons_;
}

void ElectronicWaveFunction::evaluate(const ParticlesVector<Spin> &electronsVector) {
    evaluate(electronsVector.positionsVector().asEigenVector());
}

void ElectronicWaveFunction::evaluate(const Eigen::VectorXd &electronPositionsVector) {
  assert(electronPositionsVector.rows() > 0);
  assert(electronPositionsVector.rows() % 3 == 0);
  electronPositionsVectorAsEigenVector_ = electronPositionsVector;
  auto *electronDriftCollectionArray = new double[electronPositionsVector.rows()];

  Eigen::VectorXd copy = electronPositionsVector; // TODO ugly - redesign

  amolqc_eloc(copy.data(), (int)numberOfElectrons_, &determinantProbabilityAmplitude_, &jastrowFactor_,
              electronDriftCollectionArray, &localEnergy_);

  electronDriftCollection_ = Eigen::Map<Eigen::VectorXd>( electronDriftCollectionArray,
                                                          electronPositionsVector.rows(), 1);
  delete[] electronDriftCollectionArray;
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

ParticlesVector<Spin> ElectronicWaveFunction::getElectronsVector(){
  return {PositionsVector(electronPositionsVectorAsEigenVector_), getSpinTypesVector()};
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

Eigen::VectorXd ElectronicWaveFunction::getNegativeLogarithmizedProbabilityDensityGradientCollection() {

  /*std::vector<int>frozenElectrons_={1,2,3,4,5, 10,11,12,13,14};
  for (std::vector<int>::const_iterator it = frozenElectrons_.begin(); it != frozenElectrons_.end(); ++it ){
    // index of the elctrons start at 0
    electronDriftCollection_.segment(((*it)-1)*3,3) = Eigen::VectorXd::Zero(3);
  }*/
  return -2.0 * electronDriftCollection_;
}

Eigen::VectorXd ElectronicWaveFunction::getInverseNegativeLogarithmizedProbabilityDensityGradientCollection() {
  //return getNegativeLogarithmizedProbabilityDensityGradientCollection();//.cwiseInverse();
  return (-2.0 * electronDriftCollection_).cwiseInverse();
  //return 0.5*getProbabilityAmplitude()*getProbabilityAmplitudeGradientCollection().cwiseInverse();
}

ParticlesVector<Element> ElectronicWaveFunction::getAtomsVector() const {
    return atomsVector_;
}

TypesVector<Spin > ElectronicWaveFunction::getSpinTypesVector() const {
  return spinTypesVector_;
}
