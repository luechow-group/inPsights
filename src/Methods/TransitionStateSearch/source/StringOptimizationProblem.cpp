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

#include "StringOptimizationProblem.h"
#include <iomanip>


StringOptimizationProblem::StringOptimizationProblem(long numberOfStates,
                                                     long numberOfCoords,
                                                     ElectronicWaveFunction &wf,
                                                     Eigen::MatrixXd &unitTangents //TODO make const
)
        : stepCounter_(0),
          numberOfStates_(numberOfStates),
          numberOfCoords_(numberOfCoords),
          wf_(wf),
          unitTangents_(unitTangents),
          values_(numberOfStates),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && (numberOfCoords_%3 == 0) );
    // at the beginning, no electron is at an nucleus
  for (unsigned long i = 0; i < wf_.getNumberOfElectrons()*numberOfStates_; ++i) {
    indicesOfElectronsNotAtNuclei_.push_back(i);
  }
}

double StringOptimizationProblem::value(const Eigen::VectorXd &x) {
  return stateValues(x).sum();
}

Eigen::VectorXd StringOptimizationProblem::stateValues(const Eigen::VectorXd &x) {

  Eigen::VectorXd stateValues(numberOfStates_);
  for (int i = 0; i < numberOfStates_; ++i) {
    Eigen::VectorXd xi(x.segment(i * numberOfCoords_, numberOfCoords_));

    valueCallCount_++;
    wf_.evaluate(xi);
    stateValues(i) = wf_.getNegativeLogarithmizedProbabilityDensity();
    //stateValues(i) = -wf_.getProbabilityDensity();
    //stateValues(i) = -wf_.getInverseNegativeLogarithmizedProbabilityDensity();
  }
  return stateValues;
}

void StringOptimizationProblem::fixGradient(Eigen::VectorXd &grad) {
  //std::cout << grad.transpose() << std::endl;
  for(unsigned i = 0; i < grad.size(); i++){
    if(grad[i] != grad[i]) {
      grad[i] = 0;
    }
  }
  //std::cout << grad.transpose() << std::endl;
}

void StringOptimizationProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  Eigen::VectorXd stateGrad(numberOfCoords_);
  Eigen::VectorXd unitTangent(numberOfCoords_);

  stateTypes_.resize(0);

  stateTypes_.push_back(StateGradientType::Fixed);
  for (int i = 1; i < numberOfStates_-1; ++i) stateTypes_.push_back(StateGradientType::OrthogonalToString);
  stateTypes_.push_back(StateGradientType::Fixed);

  // TODO increase efficiency by storing/updating the state values at an earlier stage
  // Identify the climbing image
  /*Eigen::VectorXd stateValuesVec = stateValues(x);
  Eigen::Index maxIdx;
  stateValuesVec.maxCoeff(&maxIdx);
  stateTypes_[maxIdx] = StateGradientType::ClimbingImage;
  */

  assert(stateTypes_.size() == numberOfStates_ );

  for (int i = 0; i < numberOfStates_; ++i) {

    gradientCallCount_++;
    wf_.evaluate(x.segment(i * numberOfCoords_, numberOfCoords_));
    values_(i) = wf_.getNegativeLogarithmizedProbabilityDensity();
    //values_(i) = -wf_.getProbabilityDensity();
    //values_(i) = -wf_.getInverseNegativeLogarithmizedProbabilityDensity();

    switch (stateTypes_[i]) {
      case StateGradientType::SimpleGradient : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
          //std::cout << i <<":"<< std::fixed << std::setw(5) << stateGrad.segment(3*8,3).transpose() << "|"<< stateGrad.segment(3*17,3).transpose() << std::endl;
        fixGradient(stateGrad);
        //stateGrad = -wf_.getProbabilityDensityGradientCollection();
        //stateGrad = -wf_.getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();
        break;
      }
      case StateGradientType::OrthogonalToString : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
          //std::cout << i <<":" << stateGrad.segment(3*0,3).transpose() << "|"<< stateGrad.segment(3*8,3).transpose() << "|"<< stateGrad.segment(3*17,3).transpose() << std::endl;
        fixGradient(stateGrad);
        //stateGrad = -wf_.getProbabilityDensityGradientCollection();
        //stateGrad = -wf_.getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();
        unitTangent = unitTangents_.row(i);
          //std::cout << i <<":" << std::fixed << std::setw(5)<< unitTangent.segment(3*0,3).transpose() << "|"<< unitTangent.segment(3*8,3).transpose() << "|"<< unitTangent.segment(3*17,3).transpose() << std::endl;
        stateGrad = stateGrad - (stateGrad.dot(unitTangent)) * unitTangent;
          //std::cout << i <<":" << std::fixed << std::setw(5) << stateGrad.segment(3*0,3).transpose() << "|"<< stateGrad.segment(3*8,3).transpose() << "|"<< stateGrad.segment(3*17,3).transpose() << std::endl;
        break;
      }
      case StateGradientType::Fixed : {
        stateGrad = Eigen::VectorXd::Zero(numberOfCoords_);
        break;
      }
      case StateGradientType::ClimbingImage : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        //stateGrad = -wf_.getProbabilityDensityGradientCollection();
        //stateGrad = -wf_.getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();

        // only if grad is small enough
        if ( stateGrad.template lpNorm<Eigen::Infinity>() < 1e-3 ) {
          // climbing
          unitTangent = unitTangents_.row(i);
          stateGrad = -stateGrad + 2 * (stateGrad.dot(unitTangent)) * unitTangent;
        } else{
          // orthogonal
          stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
          //stateGrad = -wf_.getProbabilityDensityGradientCollection();
          //stateGrad = -wf_.getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();
          unitTangent = unitTangents_.row(i);
          stateGrad = stateGrad - (stateGrad.dot(unitTangent)) * unitTangent;
        }

        break;
      }
    }
    grad.segment(i * numberOfCoords_, numberOfCoords_) = stateGrad;
  }
  //fixGradient(grad);
}


void StringOptimizationProblem::putElectronsIntoNuclei(Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  assert( x.size() == wf_.getNumberOfElectrons()*3 * numberOfStates_ && "Number of dimensions must be identical and multiple of 3");

  auto atomsVector = wf_.getAtomsVector();
  auto numberOfNuclei = atomsVector.numberOfEntities();
  auto numberOfElectrons = wf_.getNumberOfElectrons();

  // iterate over electrons that were not at nuclei in the last step
  for(auto i : indicesOfElectronsNotAtNuclei_){
    unsigned long closestNucleusIdx=0;
    double smallestDistance = std::numeric_limits<double>::infinity();

    // iterate over all nuclei and find the index of the electron with the smallest distance
    for(unsigned long j = 0; j < numberOfNuclei; ++j){
      double distance = (atomsVector[j].position()-x.segment(i*3,3)).norm();
      if(distance < smallestDistance ) {
        smallestDistance = distance;
        closestNucleusIdx = j;
      }
    }
    // check the electron with the smallest distance is close than the threshold
    double threshold = 0.05;
    if (smallestDistance <= threshold){
      //TODO PROPER?
      //Eigen::Block<Eigen::VectorXd, i*3, 0>(x.derived(), 0, 0) = atomsVector[closestNucleusIdx].position();
      x.segment(i*3,3) = atomsVector[closestNucleusIdx].position();

      //save the electron index in the indicesOfElectronsAtNuclei_ vector
      indicesOfElectronsAtNuclei_.push_back(i);
      std::sort(indicesOfElectronsAtNuclei_.begin(),indicesOfElectronsAtNuclei_.end());

      // recalulate gradient
      gradient(x,grad);
      gradientResetQ = true;
    }
  }
  // now, remove the electrons from the indicesOfElectronsNotAtNuclei_ vector
  for(auto i : indicesOfElectronsAtNuclei_) {
    indicesOfElectronsNotAtNuclei_.erase(std::remove(indicesOfElectronsNotAtNuclei_.begin(),
                                                     indicesOfElectronsNotAtNuclei_.end(), i),
                                         indicesOfElectronsNotAtNuclei_.end());
  }
}

bool StringOptimizationProblem::callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) {
  stepCounter_++;
  gradientResetQ = false;
  putElectronsIntoNuclei(x, grad); //gradientQ could be true now

  std::cout << "(" << std::setw(2) << state.iterations << ")"
            << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
            << " xDelta = " << std::setw(8) << state.xDelta
            << " gradInfNorm = " << std::setw(8) << state.gradNorm
            << std::endl;

  for (auto & it : indicesOfElectronsNotAtNuclei_) std::cout << it << " ";
  std::cout << std::endl;
  for (auto & it : indicesOfElectronsAtNuclei_) std::cout << it << " ";
  std::cout << std::endl;

  return true;
}

std::vector<unsigned long> StringOptimizationProblem::getIndicesOfElectronsNotAtNuclei() {
  return indicesOfElectronsNotAtNuclei_;
}

std::vector<unsigned long> StringOptimizationProblem::getIndicesOfElectronsAtNuclei() {
  return indicesOfElectronsAtNuclei_;
}


Eigen::VectorXd StringOptimizationProblem::getNucleiPositions() const{
  return wf_.getAtomsVector().positionsVector().asEigenVector();
}

/*
bool StringOptimizationProblem::callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) {
  gradientResetQ = false;
  putElectronsIntoNuclei(x, grad); //gradientQ could be true now

  optimizationPath_.append(ElectronsVector(x, wf_.getSpinTypesVector().asEigenVector()));

  std::cout << "(" << std::setw(2) << state.iterations << ")"
            << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
            << " xDelta = " << std::setw(8) << state.xDelta
            << " gradInfNorm = " << std::setw(8) << state.gradNorm
            << std::endl;
  std::cout << "value calls: " <<  valueCallCount_ << ", gradient calls:" << gradientCallCount_ << std::endl;

  for (auto & it : indicesOfElectronsNotAtNuclei_) std::cout << it << " ";
  std::cout << std::endl;
  for (auto & it : indicesOfElectronsAtNuclei_) std::cout << it << " ";
  std::cout << std::endl;

  return true;
}
*/
