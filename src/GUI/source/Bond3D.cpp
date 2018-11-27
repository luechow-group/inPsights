
#include <Bond3D.h>

#include "Bond3D.h"

#include "GuiHelper.h"
#include "cmath"

// bond is internally attached to the src root
Bond3D::Bond3D(const Atom3D &src, const Atom3D &dest)
  : DividedCylinder(src.parentEntity(),
                    {GuiHelper::QColorFromElementType(src.getElementType()),
                     GuiHelper::QColorFromElementType(dest.getElementType())},
                    {src.getLocation(),
                     dest.getLocation()},
                    2.4f/40.0f*std::exp(-0.1f*(src.getLocation()-dest.getLocation()).length()), 0.25f)
                    {}

Bond3D::Bond3D(Qt3DCore::QEntity *root, const Atom &src, const Atom &dest)
    : DividedCylinder(root,
                      {GuiHelper::QColorFromElementType(src.type()),
                       GuiHelper::QColorFromElementType(dest.type())},
                      {GuiHelper::toQVector3D(src.position()),
                       GuiHelper::toQVector3D(dest.position())},
                      2.4f/40.0f*std::exp(-0.1f*(
                                                        GuiHelper::toQVector3D(src.position()) -
                                      GuiHelper::toQVector3D(dest.position())).length()), 0.25f)
                              {}
