// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Cutoff.h"
#include <cmath>
#include "SOAPSettings.h"

using namespace SOAP;

bool Cutoff::withinCutoffRadiusQ(double distance) {
    return distance < Cutoff::settings.radius();
}

double Cutoff::getWeight(double distanceFromExpansionCenter) {
    const auto innerPlateauRadius = Cutoff::innerPlateauRadius();
    const auto cutoffWidth = Cutoff::settings.width();
    const auto cutoffRadius = Cutoff::settings.radius();

    if (distanceFromExpansionCenter <= innerPlateauRadius)
        return 1.;
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + cos( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth) );
    else
        return 0.;
};


double Cutoff::getWeight(const Eigen::Vector3d& position,
                         const Eigen::Vector3d& expansionCenter) {
    return getWeight(distance(position, expansionCenter));
}

Eigen::Vector3d Cutoff::getWeightGradient(const Eigen::Vector3d&position ) {
    const auto innerPlateauRadius = Cutoff::innerPlateauRadius();
    const auto cutoffWidth = Cutoff::settings.width();
    //const auto & centerWeight = ExpansionSettings::Cutoff::centerWeight;
    const auto cutoffRadius = Cutoff::settings.radius();

    double distanceFromExpansionCenter =position .norm();
    Eigen::Vector3d direction =position .normalized();

    if (distanceFromExpansionCenter <= innerPlateauRadius || distanceFromExpansionCenter > cutoffRadius)
        return Eigen::Vector3d::Zero();
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + sin( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth)*M_PI/cutoffWidth )*direction;
    else
        return Eigen::Vector3d::Zero(); //TODO is this the correct behavior?
};

double Cutoff::distance(const Eigen::Vector3d &position,
                                const Eigen::Vector3d &expansionCenter){
    return (position-expansionCenter).eval().norm();
}
