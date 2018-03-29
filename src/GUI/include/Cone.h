//
// Created by heuer on 03.08.17.
//

#ifndef AMOLQCPP_CONE_H
#define AMOLQCPP_CONE_H

#include <Qt3DExtras/QConeMesh>
#include "Abstract3dObject.h"



class Cone : public Abstract3dObject {

public:
    Cone(const Cone& cone);
    Cone(Qt3DCore::QEntity *root, QColor color,
         const std::pair<QVector3D, QVector3D>& pair,
         const float bottomRadius,
         const float topRadius = 0.0f,
         const float alpha = 1.0f);

    ~Cone() {};

    float getBottomRadius() const { return bottomRadius_; };

    void setBottomRadius(const float bottomRadius) {
      bottomRadius_ = bottomRadius;
      mesh_->setBottomRadius(bottomRadius);
    };

    void setTopRadius(const float topRadius) {
      topRadius_ = topRadius;
      mesh_->setTopRadius(topRadius);
    };

    float getLength() const { return length_; };
    QVector3D getStart() const{ return start_; };
    QVector3D getEnd() const{ return end_; };
    QVector3D getDifference() const{ return difference_; };

private:
    void rotateToOrientation(const QVector3D &orientation);

    float topRadius_, bottomRadius_, length_;
    QVector3D start_, end_, difference_;
    Qt3DExtras::QConeMesh *mesh_;
};

#endif //AMOLQCPP_CONE_H
