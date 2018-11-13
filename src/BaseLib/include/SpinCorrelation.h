//
// Created by heuer on 25.10.18.
//

#ifndef AMOLQCPP_ELECTRONPAIRANALYSIS_H
#define AMOLQCPP_ELECTRONPAIRANALYSIS_H

#include <Eigen/Core>
#include <ParticlesVector.h>

namespace SpinCorrelation{
    template<typename Type>
    Eigen::MatrixXi spinCorrelations(const TypesVector<Type> &tv){
        Eigen::MatrixXi S = Eigen::MatrixXi::Zero(tv.numberOfEntities(),tv.numberOfEntities());

        for (Eigen::Index i = 0; i < S.rows(); i++)
            for (Eigen::Index j = i + 1; j < S.cols(); j++)
                S(i,j) = tv[i]==tv[j]? 1 : -1;

        // symmetrization
        S = S.selfadjointView<Eigen::Upper>();

        return S;
    };

}

#endif //AMOLQCPP_ELECTRONPAIRANALYSIS_H
