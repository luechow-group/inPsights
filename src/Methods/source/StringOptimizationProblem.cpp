//
// Created by heuer on 13.04.17.
//

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
}

double StringOptimizationProblem::value(const Eigen::VectorXd &x) {

  /*double value = 0;
  for (int i = 0; i < numberOfStates_; ++i) {
    Eigen::VectorXd xi(x.segment(i * numberOfCoords_, numberOfCoords_));

    valueCallCount_++;
    wf_.evaluate(xi);
    value += wf_.getNegativeLogarithmizedProbabilityDensity();
  }
  return value;
   *///TODO delete
  //std::cout << "StringOptimizationProblem::value(const Eigen::VectorXd &x)" << std::endl;
  return stateValues(x).sum();
}

Eigen::VectorXd StringOptimizationProblem::stateValues(const Eigen::VectorXd &x) {

  Eigen::VectorXd stateValues(numberOfStates_);
  for (int i = 0; i < numberOfStates_; ++i) {
    Eigen::VectorXd xi(x.segment(i * numberOfCoords_, numberOfCoords_));

    //std::cout << "single state " << i << std::endl;
    //if (std::abs(x(0)) > 100){
    //  std::cout << "assert!" << std::endl;
    //}
    valueCallCount_++;
    wf_.evaluate(xi);
    stateValues(i) = wf_.getNegativeLogarithmizedProbabilityDensity();
    //std::cout << " value "<< stateValues(i) << std::endl;
  }
  return stateValues;
}

void StringOptimizationProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  Eigen::VectorXd stateGrad(numberOfCoords_);
  Eigen::VectorXd unitTangent(numberOfCoords_);


  stateTypes_.resize(0);

  stateTypes_.push_back(StateGradientType::SimpleGradient);
  for (int i = 1; i < numberOfStates_-1; ++i) {
    stateTypes_.push_back(StateGradientType::OrthogonalToString);
  }
  stateTypes_.push_back(StateGradientType::SimpleGradient);

  assert(stateTypes_.size() == numberOfStates_ );


  for (int i = 0; i < numberOfStates_; ++i) {

    gradientCallCount_++;
    wf_.evaluate(x.segment(i * numberOfCoords_, numberOfCoords_));
    values_(i) = wf_.getNegativeLogarithmizedProbabilityDensity();

    switch (stateTypes_[i]) {
      case StateGradientType::SimpleGradient : {
        stateGrad = stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        break;
      }
      case StateGradientType::OrthogonalToString : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        unitTangent = unitTangents_.row(i);

        stateGrad = stateGrad - (stateGrad.dot(unitTangent)) * unitTangent;
        break;
      }
      case StateGradientType::Fixed : {
        stateGrad = Eigen::VectorXd::Zero(numberOfCoords_);
        break;
      }
      case StateGradientType::ClimbingImage : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        unitTangent = unitTangents_.row(i);

        stateGrad = -stateGrad + 2 * (stateGrad.dot(unitTangent)) * unitTangent;
        break;
      }
    }
    grad.segment(i * numberOfCoords_, numberOfCoords_) = stateGrad;
  }
}

bool StringOptimizationProblem::callback(cppoptlib::Criteria<double> &state, Eigen::VectorXd &x) {
  stepCounter_++;
  std::cout << "(" << std::setw(2) << state.iterations << ")"
            << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
            << " xDelta = " << std::setw(8) << state.xDelta
            << " gradNorm = " << std::setw(8) << state.gradNorm
            << std::endl;
  return true;
}

