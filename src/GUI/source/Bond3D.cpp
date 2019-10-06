
#include <Bond3D.h>
#include "GuiHelper.h"
#include "cmath"

// bond is internally attached to the src root
Bond3D::Bond3D(Qt3DCore::QEntity *root, const Atom3D &src, const Atom3D &dest) //TODO what if dest gets deleted
        : DividedCylinder(root,
                          {GuiHelper::QColorFromType<Element>(src.type()),
                           GuiHelper::QColorFromType<Element>(dest.type())},
                          GuiHelper::sphericalSurfacePositionPair(
                                  src.position(), src.getRadius(),
                                  dest.position(), dest.getRadius()),
                                  0.04f,
                          0.25f) {}
