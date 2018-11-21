//
// Created by Michael Heuer on 21.11.18.
//

#include <GuiHelper.h>

QVector3D GuiHelper::midPointVector(std::pair<QVector3D, QVector3D> qVector3Dpair) {
    return qVector3Dpair.first + (qVector3Dpair.second - qVector3Dpair.first) / 2.0;
}

QColor GuiHelper::QColorFromElementType(const Element &elementType) {
        return {int(Elements::ElementInfo::color(elementType).R),
                int(Elements::ElementInfo::color(elementType).G),
                int(Elements::ElementInfo::color(elementType).B)};
}
