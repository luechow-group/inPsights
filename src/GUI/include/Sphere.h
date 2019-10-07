//
// Created by heuer on 06.12.16.
//

#ifndef INPSIGHTS_SPHERE_H
#define INPSIGHTS_SPHERE_H

#include <Qt3DExtras/QSphereMesh>
#include "Abstract3dObject.h"

class Sphere : public Abstract3dObject{
public:
    Sphere(Qt3DCore::QEntity *root, QColor color, QVector3D location, float radius, float alpha = 1.0f);

    float getRadius() const { return radius_;};

    void setRadius(const float radius) {
        radius_ = radius;
        mesh_->setRadius(radius);
    };

    friend std::ostream& operator<< (std::ostream& os, const Sphere& obj) {
        auto color = obj.color();
        auto center = obj.transform->translation();

        os << "<transform translation='"
           << center[0] << ","
           << center[1] << ","
           << center[2] << "'>\n";
        os << "<shape><appearance><material diffuseColor='"
           << color.red() << " "
           << color.green() << " "
           << color.blue()
           << "' transparency='" << obj.material->alpha() << "'></material></appearance>\n";

        os << "<sphere radius='"
           << obj.getRadius()
           << "'></sphere>\n";
        os <<"</shape></transform>\n\n";
        return os;
    }

public slots:
    void onHighlighted(bool highlightQ);
    void onSelected(bool selectedQ);

private:
    bool highlightedQ_, selectedQ_;
    float radius_;
    Qt3DExtras::QSphereMesh* mesh_;

    void update();
};

#endif //INPSIGHTS_SPHERE_H
