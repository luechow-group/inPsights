//
// Created by leonard on 22.02.18.
//

#ifndef AMOLQCPP_TESTPROBLEMS_H
#define AMOLQCPP_TESTPROBLEMS_H

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

        void hessian(const Eigen::VectorXd&x, Eigen::MatrixXd &hessian) override {
            for (int i = 0; i < x.size(); i++){
                for (int j = 0; j < x.size(); j++){
                    hessian(i,j) = 0;
                }
            }
        };

        bool callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) override{
            std::cout << "(" << std::setw(2) << state.iterations << ")"
                      << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
                      << " xDelta = " << std::setw(8) << state.xDelta
                      << " gradInfNorm = " << std::setw(8) << state.gradNorm
                      << "   " << x.transpose()
                      << std::endl;
            return true;
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

        void hessian(const Eigen::VectorXd&x, Eigen::MatrixXd &hessian) override {
            for (int i = 0; i < x.size(); i++){
                for (int j = 0; j < x.size(); j++){
                    if (i == j){
                        hessian(i,j) = 2;
                    }
                    else{
                        hessian(i,j) = 0;
                    }
                }
            }
        };

        bool callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) override{
            std::cout << "(" << std::setw(2) << state.iterations << ")"
                      << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
                      << " xDelta = " << std::setw(8) << state.xDelta
                      << " gradInfNorm = " << std::setw(8) << state.gradNorm
                      << "   " << x.transpose()
                      << std::endl;
            return true;
        };
    };
}

#endif //AMOLQCPP_TESTPROBLEMS_H
