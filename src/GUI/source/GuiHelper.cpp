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

#include <GuiHelper.h>

QVector3D GuiHelper::midPointVector(std::pair<QVector3D, QVector3D> qVector3Dpair) {
    return qVector3Dpair.first + (qVector3Dpair.second - qVector3Dpair.first) / 2.0;
}

QVector3D GuiHelper::toQVector3D(const Eigen::Vector3f &vec) {
    return {vec.x(),vec.y(), vec.z()};
}

QVector3D GuiHelper::toQVector3D(const Eigen::Vector3d &vec) {
    return toQVector3D(Eigen::Vector3f(vec.cast<float>()));
}


template<>
QColor GuiHelper::QColorFromType<Element>(const Element &type) {
    return {int(Elements::ElementInfo::color(type).R),
            int(Elements::ElementInfo::color(type).G),
            int(Elements::ElementInfo::color(type).B)};
}

template<>
QColor GuiHelper::QColorFromType<Spin>(const Spin &type) {
    switch (type) {
        case Spin::alpha:
            return Qt::red;
        case Spin::beta:
            return Qt::blue;
        case Spin::none:
            return Qt::darkMagenta;
        default:
            return Qt::darkMagenta;
    }
}

template<>
float GuiHelper::radiusFromType<Element>(const Element &type) {
    return static_cast<float>(Elements::ElementInfo::vdwRadius(type)/10.0f);
}

template<>
float GuiHelper::radiusFromType<Spin>(const Spin &type) {
    // choose electron size relative to an hydrogen atom
    return GuiHelper::radiusFromType<Element>(Element::H)/4.0f;
}

std::pair<QVector3D, QVector3D> GuiHelper::sphericalSurfacePositionPair(
        const Eigen::Vector3d& position1, double radius1,
        const Eigen::Vector3d& position2, double radius2){

    Eigen::Vector3d v12 = position2-position1;
    float distance =  v12.norm();
    Eigen::Vector3d p1 = position1 + v12 * radius1 / distance;
    Eigen::Vector3d p2 = position2 - v12 * radius2 / distance;

    return {GuiHelper::toQVector3D(p1), GuiHelper::toQVector3D(p2)};
}