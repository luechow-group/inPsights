//
// Created by heuer on 24.05.17.
//

#include "Atom.h"

Atom::Atom(Elements::ElementType elementType, Eigen::Vector3d coordinates)
        : elementType_(elementType),
          coordinates_(coordinates)
{};