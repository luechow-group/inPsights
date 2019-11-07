/* Copyright (C) 2018-2019 Michael Heuer.
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

#ifndef INPSIGHTS_RADIALGAUSSIANBASIS_H
#define INPSIGHTS_RADIALGAUSSIANBASIS_H

#include "SpecialMathFunctions/Gaussian.h"
#include <Eigen/Core>
#include <vector>

namespace SOAP {
    class RadialBasis {
    public:
        explicit RadialBasis();

        double operator()(double r, unsigned n) const;

        Eigen::MatrixXd computeCoefficients(double centerToNeighborDistance, double neighborSigma) const;

        double sigmaBasisFunction(unsigned n) { return basis_[n - 1].sigma(); };

    private:
        std::vector<Gaussian> createBasis();

        Eigen::MatrixXd Smatrix() const;

        Eigen::MatrixXd radialTransform() const;

        Eigen::MatrixXd Sab(unsigned nmax) const;

        Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);

        std::vector<double> calculateIntegrals(double ai, double ri, double rho_ik, double beta_ik) const;

        std::vector<Gaussian> basis_;
        Eigen::MatrixXd Sab_, radialTransform_;
    };
}

#endif //INPSIGHTS_RADIALGAUSSIANBASIS_H
