//
// Created by heuer on 06.12.16.
//

#ifndef TEST_SPHERE_H
#define TEST_SPHERE_H

#include <Qt3DExtras/QSphereMesh>
#include "Abstract3dObject.h"

class Sphere : public Abstract3dObject{
public:
    Sphere(Qt3DCore::QEntity *root, QColor color, QVector3D location, float radius);

    float getRadius() const { return radius_;};

    void setRadius(const float radius) {
        radius_ = radius;
        mesh_->setRadius(radius);
    };

public slots:
    void onHighlighted(bool highlightQ);
    void onSelected(bool selectedQ);

private:
    bool highlightedQ_, selectedQ_;
    float radius_;
    Qt3DExtras::QSphereMesh* mesh_;

    void update();
};

#endif //TEST_SPHERE_H
