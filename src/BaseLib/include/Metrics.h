//
// Created by Michael Heuer on 25.04.18.
//

#ifndef AMOLQCPP_METRICS_H
#define AMOLQCPP_METRICS_H


#include <Eigen/Core>
#include "PositionsVector.h"

namespace Metrics{

    template <int Norm = 2>
    double distance(const Eigen::Vector3d& position1,
            const Eigen::Vector3d& position2){
        return (position1-position2).lpNorm<Norm>();
    }

    template <int Norm = 2>
    double distance(const PositionsVector& positions1,
                    const PositionsVector& positions2){
        assert(positions1.numberOfEntities() == positions2.numberOfEntities()
               && "Both PositionVectors must have the same size.");
        return (positions1.asEigenVector()-positions2.asEigenVector()).lpNorm<Norm>();
    }

    template <int Norm = 2>
    Eigen::VectorXd positionalNormsVector(const PositionsVector &positions1, const PositionsVector &positions2) {
        assert(positions1.numberOfEntities() == positions2.numberOfEntities()
               && "Both PositionVectors must have the same size.");
        Eigen::VectorXd vec(positions1.numberOfEntities());

        for (size_t i = 0; i < positions1.numberOfEntities(); ++i) {
            vec[i] = (positions1[i]-positions2[i]).lpNorm<Norm>();
        }
        return vec;
    }

    template <int Norm = 2>
    Eigen::MatrixXd positionalDistances(const PositionsVector &positions){
        Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions.numberOfEntities(),positions.numberOfEntities());
        for (size_t i = 0; i < d.rows(); i++)
            for (size_t j = i + 1; j < d.cols(); j++)
                d(i,j) = distance<Norm>(positions[i], positions[j]);

        // symmetrization
        return d.selfadjointView<Eigen::Upper>();
    };

    template <int Norm = 2>
    Eigen::MatrixXd positionalDistances(const PositionsVector &positions1, const PositionsVector &positions2){
        assert(positions1.numberOfEntities() == positions2.numberOfEntities()
               && "Both PositionVectors must have the same size.");

        Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions1.numberOfEntities(),positions1.numberOfEntities());

        for (size_t i = 0; i < d.rows(); i++)
            for (size_t j = 0; j < d.cols(); j++)
                d(i,j) = distance<Norm>(positions1[i], positions2[j]);

        return d;
    };
}

#endif //AMOLQCPP_METRICS_H
