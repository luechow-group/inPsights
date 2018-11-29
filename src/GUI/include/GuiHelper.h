//
// Created by Michael Heuer on 25.12.16.
//

#ifndef INPSIGHTS_HELPER_H
#define INPSIGHTS_HELPER_H

#include <QColor>
#include <QVector3D>
#include <ElementInfo.h>
#include <Eigen/Core>
#include <SpinType.h>

namespace GuiHelper {
    template <typename Type>
    QColor QColorFromType(const Type &type) {
        assert(std::is_integral<Type>() && "Type must be an integer.");

        auto intType = static_cast<int>(type);

        if (intType > 0)
            return QColorFromType<Element>(Element(intType));
        else
            return QColorFromType<Spin>(Spin(intType));
    }

    template<>
    QColor QColorFromType<Spin>(const Spin &type);

    template<>
    QColor QColorFromType<Element>(const Element &type);

    QVector3D midPointVector(std::pair<QVector3D, QVector3D> qVector3Dpair);

    QVector3D toQVector3D(const Eigen::Vector3f &vec);

    QVector3D toQVector3D(const Eigen::Vector3d &vec);
}

#endif //INPSIGHTS_HELPER_H
