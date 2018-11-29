//
// Created by Michael Heuer on 21.11.18.
//

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
QColor GuiHelper::QColorFromType<Element>(const Element &type) {
    return {int(Elements::ElementInfo::color(type).R),
            int(Elements::ElementInfo::color(type).G),
            int(Elements::ElementInfo::color(type).B)};
}
