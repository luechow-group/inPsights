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

#ifndef INPSIGHTS_COULOMBPOTENTIAL_H
#define INPSIGHTS_COULOMBPOTENTIAL_H

#include <Eigen/Core>
#include "ParticlesVector.h"
#include "Metrics.h"
#include "NaturalConstants.h"
#include "EigenYamlConversion.h"

namespace CoulombPotential {

    template<typename Type>
    Eigen::MatrixXd energies(const ParticlesVector<Type> &pv, bool atomicUnits = true){
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(pv.numberOfEntities(),pv.numberOfEntities());

        for (Eigen::Index i = 0; i < V.rows()-1; i++)
            for (Eigen::Index j = i + 1; j < V.cols(); j++)
                V(i,j) = pv[i].charge() * pv[j].charge() / Metrics::distance(pv[i].position(), pv[j].position());

        // symmetrization
        V = V.selfadjointView<Eigen::Upper>();

        // optionally calculate in SI units
        if (!atomicUnits) V * std::pow(Constant::elementaryCharge,2)/(4.0*M_PI*Constant::electricConstant);

        return V;
    };

    template<typename Type1,typename Type2>
    Eigen::MatrixXd energies(const ParticlesVector<Type1> &pv1,
                             const ParticlesVector<Type2> &pv2, bool atomicUnits = true){
        Eigen::MatrixXd V(pv1.numberOfEntities(),pv2.numberOfEntities());

        for (Eigen::Index i = 0; i < V.rows(); i++)
            for (Eigen::Index j = 0; j < V.cols(); j++)
                V(i,j) = pv1[i].charge() * pv2[j].charge() / Metrics::distance(pv1[i].position(), pv2[j].position());

        // optionally calculate in SI units
        if (!atomicUnits) V * (std::pow(Constant::elementaryCharge,2) /(4.0*Constant::pi*Constant::electricConstant));

        return V;
    };
};

#endif //INPSIGHTS_COULOMBPOTENTIAL_H
