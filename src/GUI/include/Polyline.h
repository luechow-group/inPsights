//
// Created by heuer on 12.05.17.
//

#include "Cylinder.h"
#include "Cone.h"

#ifndef AMOLQCGUI_POLYLINE_H
#define AMOLQCGUI_POLYLINE_H

class Polyline : public Abstract3dObject{

public:
    Polyline(Qt3DCore::QEntity *root, QColor color, const std::vector<QVector3D> points, const float radius,
             bool arrowTipQ = false);
    ~Polyline(){};

private:
    std::vector<QVector3D> points_;
    std::vector<Cylinder*> cylinders_;
    float radius_;
    float totalArcLength_;
    bool arrowTipQ_;
    Cone* arrowTip_;
};

#endif //AMOLQCGUI_POLYLINE_H
