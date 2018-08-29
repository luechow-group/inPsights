//
// Created by Michael Heuer on 10.04.18.
//

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
    positionsVector.all().translate(-center);
};

void PositionsVectorTransformer::rotateAroundAxis(PositionsVector &p, double angle,
                                                  const Eigen::Vector3d &axisStart,
                                                  const Eigen::Vector3d &axisEnd){
    auto rotMat = rotationMatrixFromQuaternion(
            quaternionFromAngleAndAxis(angle, axisEnd - axisStart));

    p.all().translate(-axisStart);

    for (unsigned i = 0; i < p.numberOfEntities(); i++)
        p.entity(i).positionsRef() = p[i].transpose()*rotMat;

    p.all().translate(axisStart);
}

void PositionsVectorTransformer::rotateAroundAxis(PositionsVector &positionsVector, double angle,
                                                  const Eigen::Vector3d &axis){
    return rotateAroundAxis(positionsVector, angle, Eigen::Vector3d::Zero(), axis);
}

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

