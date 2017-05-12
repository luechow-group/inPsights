#include "Bond.h"

#include "Helper.h"
#include "cmath"
// bond is internally attached to the src root

Bond::Bond(const Atom &src, const Atom &dest)
  : DividedCylinder(src.parentEntity(),
                    {QColorFromElementType(src.getElementType()),
                     QColorFromElementType(dest.getElementType())},
                    {src.getLocation(),
                     dest.getLocation()},
                    2.4f/40.0f*std::exp(-0.1f*(src.getLocation()-dest.getLocation()).length())),
    src_(src),
    dest_(dest) {}