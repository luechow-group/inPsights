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

#ifndef INPSIGHTS_POSITIONSVECTORTRANSFORMER_H
#define INPSIGHTS_POSITIONSVECTORTRANSFORMER_H

#include "PositionsVector.h"
#include <Eigen/Eigenvalues>
#include <vector>

namespace PositionsVectorTransformer{

    struct AngleAxis{
        AngleAxis(double angle,Eigen::Vector3d axis)
                : angle_(angle), axis_(axis){};
        double angle_;
        Eigen::Vector3d axis_;
    };

    void translateCenterOfMassToOrigin(PositionsVector& positionsVector);

    Eigen::Matrix3d rotationMatrixFromQuaternion(const Eigen::Vector4d &q);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector,
                                   const Eigen::VectorXd& weights);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector);

    Eigen::Vector4d quaternionFromAngleAndAxis(double angle, const Eigen::Vector3d &axis);

    AngleAxis quaternionToAngleAndAxis(const Eigen::Vector4d &quaternion);
};


#endif //INPSIGHTS_POSITIONSVECTORTRANSFORMER_H
