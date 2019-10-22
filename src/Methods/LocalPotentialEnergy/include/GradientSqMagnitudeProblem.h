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

#ifndef INPSIGHTS_GRADIENTSQMAGNITUDEPROBLEM_H
#define INPSIGHTS_GRADIENTSQMAGNITUDEPROBLEM_H

#include <problem.h>
#include <iomanip>
#include "ToString.h"

using namespace ToString;

template<class ProblemType>
class GradientSqMagnitudeProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    explicit GradientSqMagnitudeProblem(ProblemType &problem)
            : problem_(problem)
            {};

    double value(const Eigen::VectorXd &x) override{
        Eigen::VectorXd grad = x;
        problem_.gradient(x,grad);

        //calculating the square of the magnitude
        double sqMagnitude = 0.0;
        for (int i = 0; i < grad.size(); i ++){
            sqMagnitude += grad(i) * grad(i);
        }
        return sqMagnitude;
    };

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override{
        Eigen::MatrixXd hessian(x.size(),x.size());
        problem_.hessian(x, hessian);

        Eigen::VectorXd gradP = x;
        problem_.gradient(x,gradP);

        for (int i = 0; i < x.size(); i++){
            grad(i) = 0;
            for (int j = 0; j < x.size(); j++){
                grad(i) += 2 * hessian(j,i) * gradP(j);
            }
        }
    };

    bool callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) override{
        if (problem_.callback(state, x, grad)){}
        return true;
    };

    ProblemType const &getProblem() const{
        return problem_;
    };

    ProblemType getProblem(){
        return problem_;
    };

private:
    ProblemType problem_;
};

#endif //INPSIGHTS_GRADIENTSQMAGNITUDEPROBLEM_H
