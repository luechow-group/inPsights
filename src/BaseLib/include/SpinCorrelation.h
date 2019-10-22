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

#ifndef INPSIGHTS_ELECTRONPAIRANALYSIS_H
#define INPSIGHTS_ELECTRONPAIRANALYSIS_H

#include <Eigen/Core>
#include <ParticlesVector.h>

namespace SpinCorrelation{
    template<typename Type>
    Eigen::MatrixXi spinCorrelations(const TypesVector<Type> &tv){
        Eigen::MatrixXi S = Eigen::MatrixXi::Zero(tv.numberOfEntities(), tv.numberOfEntities());

        for (Eigen::Index i = 0; i < S.rows(); i++)
            for (Eigen::Index j = i + 1; j < S.cols(); j++)
                S(i,j) = tv[i]==tv[j]? 1 : -1;

        // symmetrization
        S = S.selfadjointView<Eigen::Upper>();

        return S;
    };

}

#endif //INPSIGHTS_ELECTRONPAIRANALYSIS_H
