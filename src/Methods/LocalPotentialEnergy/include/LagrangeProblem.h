//
// Created by leonard on 22.02.18.
//

#ifndef AMOLQCGUI_LAGRANGEPROBLEM_H
#define AMOLQCGUI_LAGRANGEPROBLEM_H

#include <problem.h>
#include <iomanip>
#include "ToString.h"

using namespace ToString;

template<class ProblemTypeProblem, class ProblemTypeConstraint>
class LagrangeProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    LagrangeProblem(ProblemTypeProblem &problem, ProblemTypeConstraint &constraint, double equality = 0)
            : equality_(equality),
              problem_(problem),
              constraint_(constraint)
            {};

    double value(const Eigen::VectorXd &x) override{
        Eigen::VectorXd x0 = x.head(x.size()-1); // vector for problem and condition
        double lambda = x[x.size()-1];

        return problem_.value(x0) + lambda*(constraint_.value(x0) - equality_);
    };

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override{
        Eigen::VectorXd x0 = x.head(x.size()-1);

        double lambda = x[x.size()-1];

        // calculating gradients of x0 elements
        Eigen::VectorXd gradP = grad.head(grad.size()-1); // gradient of problem
        Eigen::VectorXd gradC = gradP; // gradient of condition

        constraint_.gradient(x0, gradC);
        problem_.gradient(x0, gradP);


        auto grad0 = gradP + lambda * gradC;

        // adding gradients of x0 elements to grad
        for (int i = 0; i < grad0.size(); i++) {
            grad(i) = grad0(i);
        }

        // adding gradient of lambda to grad
        grad(grad.size()-1) = constraint_.value(x0) - equality_;
    };

    bool callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) override{
        Eigen::VectorXd x0 = x.head(x.size()-1);
        if (problem_.callback(state, x0, grad)){}
        std::cout << unsignedLongToString(state.iterations,4)
                  << " | Lx= " << doubleToString(value(x),2,6)
                  << " | px= " << doubleToString(problem_.value(x0),5,0)
                  << " | cx= " << doubleToString(constraint_.value(x0) - equality_,2,6)
                  << " | ld= " << doubleToString(x(x.size()-1),5,0)
                  << " | xD= " << doubleToString(state.xDelta,5,0)
                  << " | gN= " << doubleToString(state.gradNorm,0,12)
                  << " | " << vectorXdToString(x0)
                  << std::endl;
        return true;
    };

    ProblemTypeProblem const &getProblem() const{
        return problem_;
    };

    ProblemTypeProblem getProblem(){
        return problem_;
    };

    ProblemTypeConstraint const &getConstraint() const{
        return constraint_;
    };

    ProblemTypeConstraint getConstraint(){
        return constraint_;
    };

private:
    ProblemTypeConstraint constraint_;
    ProblemTypeProblem problem_;
    double equality_;
};

#endif //AMOLQCGUI_LAGRANGEPROBLEM_H
