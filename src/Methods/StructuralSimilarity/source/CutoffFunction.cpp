//
// Created by Michael Heuer on 03.05.18.
//

#include "CutoffFunction.h"
#include <cmath>
#include "ExpansionSettings.h"

bool CutoffFunction::withinCutoffRadiusQ(double distance) {
    return distance < ExpansionSettings::Cutoff::radius;
}


double CutoffFunction::getWeight(double distanceFromExpansionCenter) {
    const auto innerPlateauRadius = ExpansionSettings::Cutoff::innerPlateauRadius();
    const auto & cutoffWidth = ExpansionSettings::Cutoff::width;
    const auto & cutoffRadius = ExpansionSettings::Cutoff::radius;

    //TODO delete centerWeight and use: 'if (0 <= distanceFromExpansionCenter...' instead?
    if (0 < distanceFromExpansionCenter && distanceFromExpansionCenter <= innerPlateauRadius)
        return 1.;
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + cos( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth) );
    else
        return 0.;
};


double CutoffFunction::getWeight(const Eigen::Vector3d& position,
                         const Eigen::Vector3d& expansionCenter) {
    return getWeight(distance(position, expansionCenter));
}

Eigen::Vector3d CutoffFunction::getWeightGradient(const Eigen::Vector3d&position ) {
    const auto innerPlateauRadius = ExpansionSettings::Cutoff::innerPlateauRadius();
    const auto & cutoffWidth = ExpansionSettings::Cutoff::width;
    const auto & centerWeight = ExpansionSettings::Cutoff::centerWeight;
    const auto & cutoffRadius = ExpansionSettings::Cutoff::radius;

    double distanceFromExpansionCenter =position .norm();
    Eigen::Vector3d direction =position .normalized();

    if (distanceFromExpansionCenter <= innerPlateauRadius || distanceFromExpansionCenter > cutoffRadius)
        return Eigen::Vector3d::Zero();
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + sin( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth)*M_PI/cutoffWidth )*direction;
};

double CutoffFunction::distance(const Eigen::Vector3d &position,
                        const Eigen::Vector3d &expansionCenter){
    return (position-expansionCenter).eval().norm();
}
