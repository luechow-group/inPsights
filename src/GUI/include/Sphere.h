// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
