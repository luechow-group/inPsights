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

#ifndef INPSIGHTS_POTENTIALPROBLEM_H
#define INPSIGHTS_POTENTIALPROBLEM_H

#include <problem.h>
#include <ParticlesVector.h>

class PotentialProblem: public cppoptlib::Problem<double,Eigen::Dynamic>{
public:

    PotentialProblem();

    explicit PotentialProblem(const AtomsVector & atomsVector);

    AtomsVector getAtoms() const;

    double value(const Eigen::VectorXd & x) override;

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override;

    void hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) override;

private:
    AtomsVector nuclei_;
    double vPotNuclei_;

    void calculateVPotNuclei();
};

#endif //INPSIGHTS_POTENTIALPROBLEM_H
