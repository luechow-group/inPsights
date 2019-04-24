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

//#include <solver/timeintegrationsolver.h>
#include <solver/gradientdescentumrigarlimitedsteplength.h>
#include <solver/gradientdescentsimplesolver.h>
#include "StringMethod.h"
#include "StringOptimizationProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/bfgssolver.h"
#include "solver/timeintegrationsolver.h"
#include "PointInterpolationGenerator.h"
#include "PenalizedLeastSquaresFitWithFixedEndsGenerator.h"
#include "PenalizedLeastSquaresFitWithLooseEndsGenerator.h"


StringMethod::StringMethod(ChainOfStates initialChain)
        ://problem.addObserver(this),
        wf_(ElectronicWaveFunction::getInstance()),
        chain_(initialChain)
{
  StringOptimizationProblem problem(chain_.statesNumber(), chain_.coordinatesNumber(), wf_, unitTangents_);
  chain_.setValues(problem.stateValues(chain_.coordinatesAsVector()));
  reparametrizeString();
}

void StringMethod::optimizeString() {
    reparametrizeString();
    discretizeStringToChain();
    calculateUnitTangents();

  unsigned maxIterations = 100000;
  unsigned iterations = 0;
    do {
        performStep();
        iterations++;
    } while (status_ == cppoptlib::Status::IterationLimit && iterations < maxIterations);
}

void StringMethod::performStep() {
    minimizeOrthogonalToString();
    reparametrizeString();
    discretizeStringToChain();
    calculateUnitTangents();
}

void StringMethod::minimizeOrthogonalToString() {
    StringOptimizationProblem problem(chain_.statesNumber(), chain_.coordinatesNumber(), wf_, unitTangents_);
    cppoptlib::TimeIntegrationSolver<StringOptimizationProblem> solver;
    //cppoptlib::BfgsnsSolver<StringOptimizationProblem> solver;
    auto crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.gradNorm = 1e-2;
    crit.iterations = 10;
    solver.setStopCriteria(crit);

    Eigen::VectorXd vec = chain_.coordinatesAsVector();
    auto oldvec = vec;

    solver.minimize(problem, vec);
    status_ = solver.status();
    chain_.storeCoordinateVectorInChain(chain_.coordinatesNumber(), chain_.statesNumber(), vec);
    // is this additional value call necessary?
    chain_.setValues(problem.stateValues(vec)); //TODO redesign?

    std::cout << "--Solver status: " << status_ << std::endl;
}

void StringMethod::reparametrizeString() {

    //arrange data for bspline generation
    Eigen::MatrixXd data(1+chain_.coordinatesNumber(), chain_.statesNumber());
    data.row(0) = chain_.values();
    data.block(1,0,chain_.coordinatesNumber(),chain_.statesNumber()) = chain_.coordinates();
    BSplines::PointInterpolationGenerator generator(data.transpose(),3,true);
    //BSplines::PenalizedLeastSquaresFitWithFixedEndsGenerator generator(data.transpose(),//Transpose to account for different data layout
    //                                                                   unsigned(chain_.statesNumber()-1),
    //                                                                   3,true,0.1);
    Eigen::VectorXi excludedDimensions(1);
    excludedDimensions << 0;

    BSplines::BSpline bs = generator.generateBSpline(1);
    arcLengthParametrizedBSpline_ = BSplines::ArcLengthParametrizedBSpline(bs,excludedDimensions);

  distributeStates();
    calculateUnitTangents();
}

void StringMethod::distributeStates() {
    uValues_.resize(chain_.statesNumber()); // THIS determines the chain of state length

    uValues_.head(1)(0) = 0;
    uValues_.tail(1)(0) = 1;

    for (int i = 1; i < chain_.statesNumber()-1 ; ++i) {
        uValues_(i) = double(i) / double(chain_.statesNumber() - 1);
    }
}

void StringMethod::discretizeStringToChain() {

    Eigen::VectorXd values(uValues_.size());
    Eigen::MatrixXd coordinates(chain_.coordinatesNumber(),uValues_.size());

    // by arc length
    for (int i = 0; i < uValues_.size(); ++i) {
        Eigen::VectorXd result = arcLengthParametrizedBSpline_.evaluate(uValues_(i));
        values(i) = result(0);
        coordinates.col(i) = result.tail(chain_.coordinatesNumber());
    }

    // by value-weighted arc length
    for (int i = 0; i < uValues_.size(); ++i) {
      Eigen::VectorXd result = arcLengthParametrizedBSpline_.evaluate(uValues_(i));
      values(i) = result(0);
      coordinates.col(i) = result.tail(chain_.coordinatesNumber());
    }

    chain_ = ChainOfStates(coordinates,values);
}

void StringMethod::calculateUnitTangents() {
    unitTangents_.resize(chain_.statesNumber(),chain_.coordinatesNumber());

    for (int i = 0; i < chain_.statesNumber() ; ++i) {
        unitTangents_.row(i) = arcLengthParametrizedBSpline_.evaluate(uValues_(i),1)
                .segment(1,chain_.coordinatesNumber()).normalized(); // exclude value dimension for tangent calculation
    }
}