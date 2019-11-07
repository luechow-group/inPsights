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

#include "PositionsVectorTransformer.h"

Eigen::Matrix3d PositionsVectorTransformer::rotationMatrixFromQuaternion(const Eigen::Vector4d &q){
    Eigen::Matrix3d rotMat = Eigen::Matrix3d::Zero();

    double q11, q22, q33;
    double q01, q02, q03, q12, q13, q23;

    q11 = q(1) * q(1);
    q22 = q(2) * q(2);
    q33 = q(3) * q(3);

    q01 = q(0) * q(1);
    q02 = q(0) * q(2);
    q03 = q(0) * q(3);
    q12 = q(1) * q(2);
    q13 = q(1) * q(3);
    q23 = q(2) * q(3);

    rotMat(0, 0) = 1 - 2*(q22 + q33);
    rotMat(1, 0) = 2*(q12 + q03);
    rotMat(2, 0) = 2*(q13 - q02);

    rotMat(0, 1) = 2*(q12 - q03);
    rotMat(1, 1) = 1 - 2*(q11 + q33);
    rotMat(2, 1) = 2*(q23 + q01);

    rotMat(0, 2) = 2*(q13 + q02);
    rotMat(1, 2) = 2*(q23 - q01);
    rotMat(2, 2) = 1 - 2*(q11 + q22);
    return rotMat;
};


void PositionsVectorTransformer::translateCenterOfMassToOrigin(PositionsVector& positionsVector){
    auto center = PositionsVectorTransformer::calculateCenterOfMass(positionsVector);
    positionsVector.translate(-center);
};

Eigen::Vector3d PositionsVectorTransformer::calculateCenterOfMass(const PositionsVector& positionsVector,
                                      const Eigen::VectorXd& weights){
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (unsigned i = 0; i < positionsVector.numberOfEntities(); i++)
        center += positionsVector[i] * weights(i);

    center /= weights.sum();
    return center;
};

Eigen::Vector3d PositionsVectorTransformer::calculateCenterOfMass(const PositionsVector& positionsVector){
    auto weights = Eigen::VectorXd::Ones(positionsVector.numberOfEntities());

    return calculateCenterOfMass(positionsVector,weights);
};

Eigen::Vector4d PositionsVectorTransformer::quaternionFromAngleAndAxis(double angle, const Eigen::Vector3d &axis){
    Eigen::Vector4d quaternion;
    quaternion <<
               cos(angle/2.),
            axis.normalized()*sin(angle/2.);
    return quaternion;
};

PositionsVectorTransformer::AngleAxis PositionsVectorTransformer::quaternionToAngleAndAxis(const Eigen::Vector4d &quaternion){
    double angle = acos(quaternion(0))*2.;
    Eigen::Vector3d axis = quaternion.tail(3);

    if(sin(angle/2.) != 0) axis /= sin(angle/2.);
    return {angle,axis};
};

