//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHOD_H
#define AMOLQCGUI_STRINGMETHOD_H

#include "ElectronicWaveFunction.h"
#include "BSpline.h"
#include "ArcLengthParametrizedBSpline.h"
#include "meta.h"

class ChainOfStates{
public:

    ChainOfStates(const Eigen::MatrixXd& coordinates, const Eigen::VectorXd& values)
            : coordinates_(coordinates),
              values_(values)
    {}
    ChainOfStates(const Eigen::MatrixXd& coordinates)
            : coordinates_(coordinates),
              values_(this->statesNumber())
    {}

    long statesNumber(){return coordinates_.cols();}
    long coordinatesNumber(){return coordinates_.rows();}

    const Eigen::MatrixXd & coordinates() { return coordinates_; }
    const Eigen::VectorXd & values() { return values_; };
    //Eigen::MatrixXd   coordinatesCopy() { return coordinates_; }
    //const Eigen::MatrixXd & coordinates() { return coordinates_; }

    Eigen::VectorXd coordinatesAsVector();

    void storeCoordinateVectorInChain(long coordinatesNumber, long statesNumber, Eigen::VectorXd &coordinateVector);

    void setCoordinates(Eigen::MatrixXd coordinates){
        coordinates_ = coordinates;
    }

    void setValues(const Eigen::VectorXd & values){
      values_ = values;
    }

private:
    //Eigen::MatrixXd coordinates_;
    Eigen::MatrixXd coordinates_;//TODO should be row major
    Eigen::VectorXd values_;
};

class StringMethod {

public:
    StringMethod(ChainOfStates chain);
    void optimizeString();

    ChainOfStates & getChain(){ return chain_; }

private:
    void minimizeOrthogonalToString();
    void performStep();
    void reparametrizeString(); // specify u values
    void distributeStates();
    void discretizeStringToChain();
    void calculateUnitTangents();

    ElectronicWaveFunction wf_;
    BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline_; // use arclength parametrized spline
    //BSplines::BSpline bSpline_;
    ChainOfStates chain_;
    Eigen::MatrixXd unitTangents_;
    Eigen::VectorXd uValues_;
    cppoptlib::Status status_;
};

/*
#include <Eigen/Core>
#include <vector>
#include "problem.h"
#include "problemobserver.h"
#include "solver/isolver.h"

template<typename SolverType>
class StringMethod : public cppoptlib::ProblemObserver{
public:
    using ProblemType = typename SolverType::ProblemType;
    using Scalar = typename SolverType::Scalar;
    using TVector = typename SolverType::TVector;

    StringMethod(unsigned numberOfStates)
            : numberOfStates(numberOfStates) {
        problem.addObserver(this);
        resetString(numberOfStates);
    }

    void resetString(unsigned numberOfStates) {
        solvers_.clear();
        for (unsigned i = 0; i < numberOfStates; ++i) {
            solvers_.push_back(SolverType());
        }
        std::cout << "Number of Solvers " << solvers_.size() << std::endl;
    }

    void evaluateString(Eigen::VectorXd &x) {
        double f = problem.value(x);
        std::cout << f << std::endl;
    }

    void stepPerformed() override{
        std::cout << "step performed" << std::endl;
    };

private:
    unsigned numberOfStates;
    ProblemType problem;
    std::vector<SolverType> solvers_;
};
*/

#endif //AMOLQCGUI_STRINGMETHOD_H
