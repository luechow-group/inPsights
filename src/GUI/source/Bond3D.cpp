
#include <Bond3D.h>
#include "GuiHelper.h"
#include "cmath"

// bond is internally attached to the src root
Bond3D::Bond3D(const Atom3D &src, const Atom3D &dest) //TODO what if dest gets deleted
        : DividedCylinder(src.parentEntity(),
                          {GuiHelper::QColorFromType<Element>(src.type()),
                           GuiHelper::QColorFromType<Element>(dest.type())},
                          {GuiHelper::toQVector3D(src.position()),
                           GuiHelper::toQVector3D(dest.position())},
                          2.4f / 40.0f * std::exp(-0.1f * (src.position() - dest.position()).norm()),
                          0.25f) {}
