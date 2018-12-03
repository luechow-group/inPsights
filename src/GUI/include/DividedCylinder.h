//
// Created by Michael Heuer on 24.12.16.
//

#ifndef INPSIGHTS_DIVIDEDCYLINDER_H
#define INPSIGHTS_DIVIDEDCYLINDER_H

#include <Qt3DCore/QEntity>
#include "Cylinder.h"

class DividedCylinder : public Qt3DCore::QEntity{
public:
  DividedCylinder(Qt3DCore::QEntity *root,
                  const std::pair<QColor,QColor>& colorPair,
                  const std::pair<QVector3D, QVector3D>& locationPair,
                  float radius,
                  float alpha = 1.0f);

 ~DividedCylinder() { /*QT manages destruction*/ };

private:
  Cylinder *srcCylinder_, *destCylinder_;
};

#endif //INPSIGHTS_DIVIDEDCYLINDER_H
