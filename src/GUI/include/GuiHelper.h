//
// Created by Michael Heuer on 25.12.16.
//

#ifndef INPSIGHTS_HELPER_H
#define INPSIGHTS_HELPER_H

#include <QColor>
#include <QVector3D>
#include <ElementInfo.h>

namespace GuiHelper {
    QColor QColorFromElementType(const Element &elementType);

    QVector3D midPointVector(std::pair<QVector3D, QVector3D> qVector3Dpair);
}

#endif //INPSIGHTS_HELPER_H
