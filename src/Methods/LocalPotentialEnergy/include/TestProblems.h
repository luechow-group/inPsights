//
// Created by leonard on 22.02.18.
//

#ifndef AMOLQCGUI_TESTPROBLEMS_H
#define AMOLQCGUI_TESTPROBLEMS_H

#include <problem.h>

namespace TestProblems{
    class testProblem: public cppoptlib::Problem<double,Eigen::Dynamic>{
    public:
        double value(const Eigen::VectorXd & x) override {
            double value = 0.0;
            for (int i = 0; i < x.size(); i++){
                value += x[i];
            }
            return value;
        };

        void gradient(const Eigen::VectorXd & x, Eigen::VectorXd &grad) override {
            for (int i = 0; i < x.size(); i++){
                grad[i] = 1;
            }
        };
    };

    class testConstraint : public cppoptlib::Problem<double,Eigen::Dynamic>{
    public:
        double value(const Eigen::VectorXd & x) override {
            double value = 0.0;
            for (int i = 0; i < x.size(); i++){
                value += x[i]*x[i];
            }
            return value;
        };

        void gradient(const Eigen::VectorXd & x, Eigen::VectorXd &grad) override {
            for (int i = 0; i < x.size(); i++){
                grad[i] = 2*x[i];
            }
        };
    };
}

#endif //AMOLQCGUI_TESTPROBLEMS_H
