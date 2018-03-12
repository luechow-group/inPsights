//
// Created by leonard on 22.02.18.
//

#ifndef AMOLQCGUI_LAGRANGEPROBLEM_H
#define AMOLQCGUI_LAGRANGEPROBLEM_H

#include <problem.h>
#include <iomanip>

template<class ProblemTypeProblem, class ProblemTypeConstraint>
class LagrangeProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    LagrangeProblem(ProblemTypeProblem &problem, ProblemTypeConstraint &constraint, double equality = 0)
            : equality_(equality),
              problem_(problem),
              constraint_(constraint)
            {};

    double value(const Eigen::VectorXd &x) override{
        auto x0 = x.block(0,0,x.size()-1,1); // vector for problem and condition
        double lambda = x[x.size()-1];

        return problem_.value(x0) + lambda*(constraint_.value(x0) - equality_);
    };

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override{
        Eigen::VectorXd x0 = x.block(0,0,x.size()-1,1);

        double lambda = x[x.size()-1];

        // calculating gradients of x0 elements
        Eigen::VectorXd gradP = grad.block(0,0,x.size()-1,1); // gradient of problem
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
        std::cout << "(" << std::setw(2) << state.iterations << ")"
                  << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
                  << " xDelta = " << std::setw(8) << state.xDelta
                  << " gradInfNorm = " << std::setw(8) << state.gradNorm
                  << "   " << x.transpose()
                  << std::endl;
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
