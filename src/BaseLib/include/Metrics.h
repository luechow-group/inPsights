//
// Created by Michael Heuer on 25.04.18.
//

#ifndef INPSIGHTS_METRICS_H
#define INPSIGHTS_METRICS_H


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
    double minimalDistance(const Eigen::Vector3d &position,
                           const std::vector<Eigen::Vector3d> &positions){
        double minDistance = -1.0;
        double distance;
        for (auto position2 : positions){
            distance = (position-position2).lpNorm<Norm>();
            if (distance < minDistance or minDistance < 0){
                minDistance = distance;
            }
        }
        return minDistance;
    }

    template <int Norm = 2>
    double averageDistance(const Eigen::Vector3d &position,
                           const std::vector<Eigen::Vector3d> &positions){
        double avgDistance = 0.0;
        for (auto position2 : positions){
            avgDistance += (position-position2).lpNorm<Norm>();
        }
        return avgDistance/positions.size();
    }

    template <int Norm = 2>
    Eigen::VectorXd positionalNormsVector(const PositionsVector &positions1, const PositionsVector &positions2) {
        assert(positions1.numberOfEntities() == positions2.numberOfEntities()
               && "Both PositionVectors must have the same size.");
        Eigen::VectorXd vec(positions1.numberOfEntities());

        for (Eigen::Index  i = 0; i < positions1.numberOfEntities(); ++i) {
            vec[i] = (positions1[i]-positions2[i]).lpNorm<Norm>();
        }
        return vec;
    }

    template <int Norm = 2>
    Eigen::MatrixXd positionalDistances(const PositionsVector &positions){
        Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions.numberOfEntities(),positions.numberOfEntities());
        for (Eigen::Index  i = 0; i < d.rows(); i++)
            for (Eigen::Index  j = i + 1; j < d.cols(); j++)
                d(i,j) = distance<Norm>(positions[i], positions[j]);

        // symmetrization
        return d.selfadjointView<Eigen::Upper>();
    };

    template <int Norm = 2>
    Eigen::MatrixXd positionalDistances(const PositionsVector &positions1, const PositionsVector &positions2){
        Eigen::MatrixXd d = Eigen::MatrixXd::Zero(positions1.numberOfEntities(),positions2.numberOfEntities());

        for (Eigen::Index  i = 0; i < d.rows(); i++)
            for (Eigen::Index j = 0; j < d.cols(); j++)
                d(i,j) = distance<Norm>(positions1[i], positions2[j]);

        return d;
    };

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    double positionalNormsVectorNorm(
            const PositionsVector &p1,
            const PositionsVector &p2) {
        assert(p1.numberOfEntities() == p2.numberOfEntities());

        Eigen::VectorXd res = Metrics::positionalNormsVector<positionalNorm>(p1, p2);
        return res.lpNorm<overallNorm>();
    }
}

#endif //INPSIGHTS_METRICS_H
