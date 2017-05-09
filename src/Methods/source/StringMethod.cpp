//
// Created by heuer on 21.04.17.
//

#include "StringMethod.h"
#include "StringOptimizationProblem.h"
#include "solver/bfgsnssolver.h"
#include "PointInterpolationGenerator.h"
#include "PenalizedLeastSquaresFitWithFixedEndsGenerator.h"
#include "PenalizedLeastSquaresFitWithLooseEndsGenerator.h"


Eigen::VectorXd ChainOfStates::coordinatesAsVector() {

  //TODO use map for efficiency
  Eigen::VectorXd vec(statesNumber()*coordinatesNumber());

  for (int i = 0; i < statesNumber(); ++i) {
    vec.segment(i*coordinatesNumber(),coordinatesNumber()) = coordinates().row(i);
  }
  //TODO use map for efficiency
  //std::cout << coordinates() << std::endl;
  //std::cout << vec.transpose() << std::endl;
  //return Eigen::Map<Eigen::VectorXd>(coordinates_.data(), statesNumber()*coordinatesNumber());
  return vec;
};


void ChainOfStates::storeVectorInChain(long statesNumber, long coordinatesNumber, Eigen::VectorXd &vec) {
  for (int i = 0; i < statesNumber; ++i) {
    coordinates().row(i) = vec.segment(i*coordinatesNumber,coordinatesNumber);
  }
  //Eigen::Map< Eigen::MatrixXd>temp(vec.data(),statesNumber,coordinatesNumber);
    //coordinates_.Map(vec.data(),statesNumber,coordinatesNumber);
    //coordinates_ = temp;
}

StringMethod::StringMethod(ChainOfStates initialChain)
        : chain_(initialChain)
{}

void StringMethod::optimizeString() {
    reparametrizeString();
    discretizeStringToChain();
    calculateUnitTangents();
    do {
        performStep();
    } while (status_ == cppoptlib::Status::IterationLimit);
}

void StringMethod::performStep() {
    minimizeOrthogonalToString();
    reparametrizeString();
    discretizeStringToChain();
    calculateUnitTangents();
    // for discretized states //it inside discretize chain -> it makes no sense to call it to a different point in time
}

void StringMethod::minimizeOrthogonalToString() {
    StringOptimizationProblem problem(chain_.statesNumber(), chain_.coordinatesNumber(), wf_, unitTangents_);
    cppoptlib::BfgsnsSolver<StringOptimizationProblem> solver;
    auto crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.gradNorm = 0.00001;
    crit.iterations = 10;
    solver.setStopCriteria(crit);

    Eigen::VectorXd vec = chain_.coordinatesAsVector();
    auto oldvec = vec;

    solver.minimize(problem, vec);
    status_ = solver.status();
    chain_.storeVectorInChain(chain_.statesNumber(), chain_.coordinatesNumber(), vec);

    //std::cout << " old "<< oldvec.transpose() << std::endl;
    //std::cout << "diff " << (oldvec-vec).transpose() << std::endl;
    //std::cout << " new " << vec.transpose() << std::endl;
    std::cout << chain_.coordinates() << std::endl;

    std::cout << "--Solver status: " << status_ << std::endl;
}

void StringMethod::reparametrizeString() {

    //TODO also fit energies, employ energy weighting
    //BSplines::PointInterpolationGenerator generator(chain_.coordinates(),3,true);
    BSplines::PenalizedLeastSquaresFitWithFixedEndsGenerator generator(chain_.coordinates(),
                                                                       unsigned(chain_.statesNumber()-1),
                                                                       3,true,0.4);
    Eigen::VectorXi excludedDimensions(0);
    //excludedDimensions << 1;

    BSplines::BSpline bs = generator.generateBSpline(1);

    arcLengthParametrizedBSpline_ = BSplines::ArcLengthParametrizedBSpline(bs,excludedDimensions);

    calculateParameterValues();
    calculateUnitTangents();
}

void StringMethod::calculateParameterValues() {
    uValues_.resize(chain_.statesNumber());

    uValues_.head(1)(0) = 0;
    uValues_.tail(1)(0) = 1;

    for (int i = 1; i < chain_.statesNumber()-1 ; ++i) {
        uValues_(i) = double(i) / double(chain_.statesNumber() - 1);
    }
}

void StringMethod::discretizeStringToChain() {

    for (int i = 0; i < chain_.statesNumber() ; ++i) {
      //std::cout << uValues_(i) << ": "
      //          << arcLengthParametrizedBSpline_.evaluate(uValues_(i)).transpose() << std::endl;

      chain_.coordinates().row(i) = arcLengthParametrizedBSpline_.evaluate(uValues_(i));
    }
}

void StringMethod::calculateUnitTangents() {
    unitTangents_.resize(chain_.statesNumber(),chain_.coordinatesNumber());

    for (int i = 0; i < chain_.statesNumber() ; ++i) {
        unitTangents_.row(i) = arcLengthParametrizedBSpline_.evaluate(uValues_(i),1).normalized();
    }
}