//
// Created by Michael Heuer on 25.12.16.
//

#ifndef AMOLQCPP_HELPER_H
#define AMOLQCPP_HELPER_H

#include <QColor>
#include <QVector3D>

#include "ElementInfo.h"


static QColor QColorFromElementType(const Element& elementType){
  return QColor(Elements::ElementInfo::color(elementType).R,
                Elements::ElementInfo::color(elementType).G,
                Elements::ElementInfo::color(elementType).B);
}

static QVector3D MidPointVector(std::pair<QVector3D,QVector3D> qVector3Dpair){
  return qVector3Dpair.first + (qVector3Dpair.second - qVector3Dpair.first) / 2.0;
}




#endif //AMOLQCPP_HELPER_H
