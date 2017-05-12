//
// Created by heuer on 12.05.17.
//

#include "Polyline.h"
#include "Helper.h"

Polyline::Polyline(Qt3DCore::QEntity *root,
                   QColor color,
                   const std::vector<QVector3D> points,
                   const float radius)
        : Abstract3dObject(root,color,MidPointVector(std::make_pair(*points.begin(),*points.end()))),
          points_(points),
          cylinders_(0),
          radius_(radius)
{

  //for (std::vector<QVector3D>::const_iterator point = points.begin()+1; point != points.end(); point++){
  for (int i = 1; i < points.size(); ++i) {
    Cylinder c (root,color,std::make_pair(points[i-1],points[i]),radius);
    cylinders_.push_back(&c);
  }

}


