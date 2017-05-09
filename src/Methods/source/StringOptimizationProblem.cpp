//
// Created by heuer on 13.04.17.
//

#include "StringOptimizationProblem.h"
#include <iostream>
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
        //chain_(chain),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && (numberOfCoords_%3 == 0) );
}

double StringOptimizationProblem::value(const Eigen::VectorXd &x) {

  double value = 0;
  for (int i = 0; i < numberOfStates_; ++i) {
    Eigen::VectorXd xi(x.segment(i * numberOfCoords_, numberOfCoords_));

    valueCallCount_++;
    wf_.evaluate(xi);

    value += wf_.getNegativeLogarithmizedProbabilityDensity();
    //chain_(i,0) = wf_.getNegativeLogarithmizedProbabilityDensity();
  }
  return value;//chain_.col(0).sum();
}

void StringOptimizationProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  Eigen::VectorXd stateGrad(numberOfCoords_);
  Eigen::VectorXd unitTangent(numberOfCoords_);


  stateTypes_.resize(0);

  stateTypes_.push_back(StateType::Simple);
  for (int i = 1; i < numberOfStates_-1; ++i) {
    stateTypes_.push_back(StateType::Orthogonal);
  }
  stateTypes_.push_back(StateType::Simple);

  assert(stateTypes_.size() == numberOfStates_ );


  for (int i = 0; i < numberOfStates_; ++i) {

    gradientCallCount_++;
    wf_.evaluate(x.segment(i * numberOfCoords_, numberOfCoords_));
    //std::cout << "geom " << i << " " << x.segment(i * numberOfCoords_, numberOfCoords_).transpose() << std::endl;
    switch (stateTypes_[i]) {
      case StateType::Simple : {
        stateGrad = stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        break;
      }
      case StateType::Orthogonal : {
        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
        unitTangent = unitTangents_.row(i);

        stateGrad = stateGrad - (stateGrad.dot(unitTangent)) * unitTangent;
        break;
      }
      case StateType::Fixed : {
        stateGrad = Eigen::VectorXd::Zero(numberOfCoords_);
        break;
      }
      case StateType::ClimbingImage : {
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
            << std::endl;
  return true;
}

