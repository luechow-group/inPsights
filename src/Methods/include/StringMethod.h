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
    ChainOfStates(Eigen::MatrixXd chain)
            : coordinates_(chain),
              values_(this->statesNumber())
    {}

    //ChainOfStates(Eigen::MatrixXd chain, std::vector<unsigned> valuePositins);

    long statesNumber(){return coordinates_.rows();}
    long coordinatesNumber(){return coordinates_.cols();}

    Eigen::MatrixXd & coordinates() { return coordinates_; }
    //Eigen::MatrixXd   coordinatesCopy() { return coordinates_; }
    //const Eigen::MatrixXd & coordinates() { return coordinates_; }

    Eigen::VectorXd coordinatesAsVector() {
        return Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(coordinates_.data(),
                                                    statesNumber()*coordinatesNumber()));
    };

    void storeVectorInChain(long statesNumber,
                            long coordinatesNumber,
                            const Eigen::VectorXd &vec//, std::vector<unsigned> valuePositions
    );

    void setCoordinates(Eigen::MatrixXd coordinates){
        coordinates_ = coordinates;
    }

private:
    Eigen::MatrixXd coordinates_;
    Eigen::VectorXd values_;
};

class StringMethod {

public:
    StringMethod(ChainOfStates chain);
    void optimizeString();

    ChainOfStates & getChain(){ return chain_; }

private:
    ElectronicWaveFunction wf_;
    //BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline_; // use arclength parametrized spline
    BSplines::BSpline bSpline_;
    ChainOfStates chain_;
    Eigen::VectorXd unitTangent_, uValues_;
    cppoptlib::Status status_;

//TODO make
private:
// public:
    void minimizeOrthogonalToString();
    void performStep();
    void reparametrizeString(); // specify u values
    void calculateParameterValues();
    void discretizeStringToChain();
    void calculateUnitTangent();
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
