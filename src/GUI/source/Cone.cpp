
#include <cmath>
#include "Cone.h"
#include "Helper.h"

Cone::Cone(Qt3DCore::QEntity *root,
           QColor color,
           const std::pair<QVector3D, QVector3D>& pair,
           const float bottomRadius,
           const float topRadius,
           const float alpha)
        : Abstract3dObject(root, QColor(), MidPointVector(pair)),
          topRadius_(topRadius),
          bottomRadius_(bottomRadius),
          start_(pair.first),
          end_(pair.second) {

  difference_ = end_ - start_;
  length_ = difference_.length();

  mesh_ = new Qt3DExtras::QConeMesh;
  mesh_->setTopRadius(topRadius);
  mesh_->hasTopEndcapChanged(true);
  mesh_->setBottomRadius(bottomRadius);
  mesh_->hasBottomEndcapChanged(true);
  mesh_->setLength(length_);
  mesh_->setRings(100);
  mesh_->setSlices(10);

  material->setAlpha(alpha);

  rotateToOrientation(difference_);

  material->setAmbient(color);

  entity->addComponent(mesh_);
}

Cone::Cone(const Cone &cone)
        : bottomRadius_(cone.getBottomRadius()),
          length_(cone.getLength()),
          start_(cone.getStart()),
          end_(cone.getEnd()),
          difference_(cone.difference_)
{
}

void Cone::rotateToOrientation(const QVector3D &orientation) {

  auto origVec = QVector3D(0, 1, 0);

  auto perpendicular = QVector3D::crossProduct(origVec, orientation);
  auto l = perpendicular.length();
  float epsilon = 0.00001;

  if (l > epsilon) {
    QVector3D axis = perpendicular.normalized();
    float angle = float(atan2f(l, QVector3D::dotProduct(origVec, orientation)) / (M_PI * 2) * 360.0f);
    transform->setRotation(QQuaternion::fromAxisAndAngle(axis, angle));
  } else if (QVector3D::dotProduct(origVec, orientation) > 0.0f) {
    // Nearly positively aligned; skip rotation, or compute
    // axis and angle using other means
  } else {
    // Nearly negatively aligned; axis is any vector perpendicular
    // to either vector, and angle is 180 degrees
    auto tmp = start_;
    start_ = end_;
    end_ = tmp;
  }
}
