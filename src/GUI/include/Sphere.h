/* Copyright (C) 2016-2019 Michael Heuer.
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

#ifndef INPSIGHTS_SPHERE_H
#define INPSIGHTS_SPHERE_H

#include <Qt3DExtras/QSphereMesh>
#include "Abstract3dObject.h"

class Sphere : public Abstract3dObject{
public:
    Sphere(Qt3DCore::QEntity *root, QColor color, QVector3D location, float radius, float alpha = 1.0f);

    float getRadius() const;

    void setRadius(const float radius);

    void addToXml (std::ostream& os, unsigned sortKey = 1) const;

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
