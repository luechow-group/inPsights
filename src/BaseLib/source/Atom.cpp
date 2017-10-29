//
// Created by heuer on 24.05.17.
//

#include "Atom.h"

Atom::Atom(Eigen::Vector3d position, Elements::ElementType elementType)
        : Particle(position),
          elementType_(elementType)
{};