//
// Created by Michael Heuer on 25.12.16.
//

#ifndef INPSIGHTS_HELPER_H
#define INPSIGHTS_HELPER_H

#include <QColor>
#include <QVector3D>

#include "ElementInfo.h"


static QColor QColorFromElementType(const Element& elementType){
  return {int(Elements::ElementInfo::color(elementType).R),
          int(Elements::ElementInfo::color(elementType).G),
          int(Elements::ElementInfo::color(elementType).B)};
}

static QVector3D MidPointVector(std::pair<QVector3D,QVector3D> qVector3Dpair){
  return qVector3Dpair.first + (qVector3Dpair.second - qVector3Dpair.first) / 2.0;
}

#endif //INPSIGHTS_HELPER_H
