/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_CONE_H
#define INPSIGHTS_CONE_H

#include <Qt3DExtras/QConeMesh>
#include "Abstract3dObject.h"

class Cone : public Abstract3dObject {

public:
    Cone(Qt3DCore::QEntity *root, QColor color,
         const std::pair<QVector3D, QVector3D>& pair,
         float bottomRadius,
         float topRadius = 0.0f,
         float alpha = 1.0f);

    float getBottomRadius() const { return bottomRadius_; };

    void setBottomRadius(const float bottomRadius) {
      bottomRadius_ = bottomRadius;
      mesh_->setBottomRadius(bottomRadius);
    };

    void setTopRadius(const float topRadius) {
      topRadius_ = topRadius;
      mesh_->setTopRadius(topRadius);
    };

    float length() const { return difference().length(); };
    QVector3D start() const{ return start_; };
    QVector3D end() const{ return end_; };
    QVector3D difference() const{ return end_-start_; };

private:
    void rotateToOrientation(const QVector3D &orientation);

    float topRadius_, bottomRadius_;
    QVector3D start_, end_;
    Qt3DExtras::QConeMesh *mesh_;
};

#endif //INPSIGHTS_CONE_H
