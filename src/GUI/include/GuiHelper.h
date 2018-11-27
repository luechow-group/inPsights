//
// Created by Michael Heuer on 25.12.16.
//

#ifndef INPSIGHTS_HELPER_H
#define INPSIGHTS_HELPER_H

#include <QColor>
#include <QVector3D>
#include <ElementInfo.h>
#include <Eigen/Core>

namespace GuiHelper {
    QColor QColorFromElementType(const Element &elementType);

    QVector3D midPointVector(std::pair<QVector3D, QVector3D> qVector3Dpair);

    QVector3D toQVector3D(const Eigen::Vector3f &vec);

    QVector3D toQVector3D(const Eigen::Vector3d &vec);
}

#endif //INPSIGHTS_HELPER_H
