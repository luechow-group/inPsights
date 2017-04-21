//
// Created by heuer on 21.04.17.
//

#include "StringMethod.h"
#include "StringOptimizationProblem.h"
#include "solver/bfgsnssolver.h"
#include "PointInterpolationGenerator.h"


void ChainOfStates::storeVectorInChain(long statesNumber, long coordinatesNumber, const Eigen::VectorXd &vec) {
    Eigen::Map<const Eigen::MatrixXd>temp(vec.data(),statesNumber,coordinatesNumber);
    //coordinates_.Map(vec.data(),statesNumber,coordinatesNumber);
    coordinates_ = temp;
}

StringMethod::StringMethod(ChainOfStates initialChain)
        : chain_(initialChain)
{}

void StringMethod::optimizeString() {
    reparametrizeString();
    discretizeStringToChain();
    calculateUnitTangent();
    do {
        performStep();
    } while (status_ == cppoptlib::Status::IterationLimit);
}

void StringMethod::performStep() {
    minimizeOrthogonalToString();
    reparametrizeString();
    discretizeStringToChain();
    // for discretized states //it inside discretize chain -> it makes no sense to call it to a different point in time
}

void StringMethod::minimizeOrthogonalToString() {
    StringOptimizationProblem problem(chain_.statesNumber(), chain_.coordinatesNumber(), wf_, unitTangent_);
    cppoptlib::BfgsnsSolver<StringOptimizationProblem> solver;
    auto crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.iterations = 20;
    solver.setStopCriteria(crit);

    Eigen::VectorXd vec = chain_.coordinatesAsVector();
    solver.minimize(problem, vec);
    status_ = solver.status();
    chain_.storeVectorInChain(chain_.statesNumber(), chain_.coordinatesNumber(), vec);

    std::cout << "--Solver status: " << status_ << std::endl;
    std::cout << chain_.coordinates() << std::endl;


}

void StringMethod::reparametrizeString() {
    // include probability density in the interpolation
    // construct (value weighted) arclength parametrized spline
    BSplines::PointInterpolationGenerator interpolator(chain_.coordinates(),3,true);
    /*Eigen::VectorXi excludedDimensions(1);
    excludedDimensions << 1;
    arcLengthParametrizedBSpline_ = BSplines::ArcLengthParametrizedBSpline(interpolator.generateBSpline(1),
                                                                           excludedDimensions);*/
    bSpline_ = interpolator.generateBSpline(1);

    calculateParameterValues();
    calculateUnitTangent();
}

void StringMethod::calculateParameterValues() {
    uValues_.resize(chain_.statesNumber());
    for (int i = 0; i < chain_.statesNumber() ; ++i) {
        uValues_(i) = double(i) / double(chain_.statesNumber() - 1);
    }
}

void StringMethod::discretizeStringToChain() {
    for (int i = 0; i < chain_.statesNumber() ; ++i) {
        chain_.coordinates().row(i) = bSpline_.evaluate(uValues_(i));
    }
    //std::cout << "spline \n"<< chain_.coordinates() << std::endl;
}

void StringMethod::calculateUnitTangent() {
    unitTangent_.resize(chain_.statesNumber()*chain_.coordinatesNumber());

    for (int i = 0; i < chain_.statesNumber() ; ++i) {
        unitTangent_.segment(i*chain_.coordinatesNumber(), chain_.coordinatesNumber()) = bSpline_.evaluate(uValues_(i),1);
    }
    unitTangent_.normalize();
}