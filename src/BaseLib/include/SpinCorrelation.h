// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ELECTRONPAIRANALYSIS_H
#define INPSIGHTS_ELECTRONPAIRANALYSIS_H

#include <Eigen/Core>
#include <ParticlesVector.h>

namespace SpinCorrelation{
    template<typename Type>
    Eigen::MatrixXi spinCorrelations(const TypesVector<Type> &tv){
        Eigen::MatrixXi S = Eigen::MatrixXi::Zero(tv.numberOfEntities(), tv.numberOfEntities());

        for (Eigen::Index i = 0; i < S.rows()-1; i++)
            for (Eigen::Index j = i + 1; j < S.cols(); j++)
                S(i,j) = tv[i]==tv[j]? 1 : -1;

        // symmetrization
        S = S.selfadjointView<Eigen::Upper>();

        return S;
    };

}

#endif //INPSIGHTS_ELECTRONPAIRANALYSIS_H
