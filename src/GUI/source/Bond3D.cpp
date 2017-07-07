#include "Bond3D.h"

#include "Helper.h"
#include "cmath"
// bond is internally attached to the src root

Bond3D::Bond3D(const Atom3D &src, const Atom3D &dest)
  : DividedCylinder(src.parentEntity(),
                    {QColorFromElementType(src.getElementType()),
                     QColorFromElementType(dest.getElementType())},
                    {src.getLocation(),
                     dest.getLocation()},
                    2.4f/40.0f*std::exp(-0.1f*(src.getLocation()-dest.getLocation()).length()), 0.75f),
    src_(src),
    dest_(dest) {}