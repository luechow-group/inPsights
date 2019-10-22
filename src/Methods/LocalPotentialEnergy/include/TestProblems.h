/* Copyright (C) 2018 Leonard Reuter.
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

#ifndef INPSIGHTS_TESTPROBLEMS_H
#define INPSIGHTS_TESTPROBLEMS_H

#include <problem.h>

namespace TestProblems{
    class TestProblem: public cppoptlib::Problem<double,Eigen::Dynamic>{
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
            /*std::cout << "(" << std::setw(2) << state.iterations << ")"
                      << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
                      << " xDelta = " << std::setw(8) << state.xDelta
                      << " gradInfNorm = " << std::setw(8) << state.gradNorm
                      << "   " << x.transpose()
                      << std::endl;*/
            return true;
        };
    };

    class TestConstraint : public cppoptlib::Problem<double,Eigen::Dynamic>{
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
            /*std::cout << "(" << std::setw(2) << state.iterations << ")"
                      << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
                      << " xDelta = " << std::setw(8) << state.xDelta
                      << " gradInfNorm = " << std::setw(8) << state.gradNorm
                      << "   " << x.transpose()
                      << std::endl;*/
            return true;
        };
    };
}

#endif //INPSIGHTS_TESTPROBLEMS_H
