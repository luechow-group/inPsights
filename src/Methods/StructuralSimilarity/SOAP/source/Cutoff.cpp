// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Cutoff.h"
#include <cmath>
#include "SOAPSettings.h"

using namespace SOAP;

bool Cutoff::withinCutoffRadiusQ(double distance) {
    return distance < Cutoff::settings.radius();
}

double Cutoff::value(double distanceFromExpansionCenter) {
    const auto innerPlateauRadius = Cutoff::innerPlateauRadius();
    const auto cutoffWidth = Cutoff::settings.width();
    const auto cutoffRadius = Cutoff::settings.radius();

    if (distanceFromExpansionCenter <= innerPlateauRadius)
        return 1.;
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + std::cos( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth) );
    else
        return 0.;
};


double Cutoff::value(const Eigen::Vector3d& position,
                     const Eigen::Vector3d& expansionCenter) {
    return value((position-expansionCenter).norm());
}

Eigen::Vector3d Cutoff::gradient(const Eigen::Vector3d&position ) {
    const auto innerPlateauRadius = Cutoff::innerPlateauRadius();
    const auto cutoffWidth = Cutoff::settings.width();
    const auto cutoffRadius = Cutoff::settings.radius();

    double distanceFromExpansionCenter = position.norm();
    Eigen::Vector3d direction = position.normalized();

    if (distanceFromExpansionCenter <= innerPlateauRadius || distanceFromExpansionCenter > cutoffRadius)
        return Eigen::Vector3d::Zero();
    else if (innerPlateauRadius < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius)
        return 0.5*( 1 + std::sin( M_PI*(distanceFromExpansionCenter-innerPlateauRadius)/cutoffWidth)*M_PI/cutoffWidth )*direction;
    else
        return Eigen::Vector3d::Zero(); //TODO is this the correct behavior?
};

Eigen::Vector3d Cutoff::gradient(const Eigen::Vector3d& position,
                                 const Eigen::Vector3d& expansionCenter) {
    return gradient(position-expansionCenter);
}
