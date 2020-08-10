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

#ifndef INPSIGHTS_METRICS_H
#define INPSIGHTS_METRICS_H


#include <Eigen/Core>
#include "PositionsVector.h"

namespace Metrics{

    template <int Norm = 2>
    double distance(const Eigen::Vector3d& position1, const Eigen::Vector3d& position2){
        return (position1-position2).lpNorm<Norm>();
    }

    template <int Norm = 2>
    double distance(const PositionsVector& positions1,
                    const PositionsVector& positions2){
        assert(positions1.numberOfEntities() == positions2.numberOfEntities()
               && "Both PositionVectors must have the same size.");
        return (positions1.asEigenVector()-positions2.asEigenVector()).lpNorm<Norm>();
    }

    Eigen::Vector3d averagedPosition(const PositionsVector& positions);

    template <int Norm = 2>
    double minimalDistance(const Eigen::Vector3d &position,
                           const std::vector<Eigen::Vector3d> &positions){
        double minDistance = (position-positions.front()).lpNorm<Norm>();
        double distance;

        for(auto positionIt = std::next(positions.begin()); positionIt != positions.end(); ++positionIt) {
            distance = (position-*positionIt).lpNorm<Norm>();
            if(distance < minDistance) {
                minDistance = distance;
            }
        }
        return minDistance;
    }

    template <int Norm = 2>
    double averageDistance(const Eigen::Vector3d &position,
                           const std::vector<Eigen::Vector3d> &positions){
        double avgDistance = 0.0;
        for (const auto& position2 : positions){
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
        for (Eigen::Index  i = 0; i < d.rows()-1; i++)
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
