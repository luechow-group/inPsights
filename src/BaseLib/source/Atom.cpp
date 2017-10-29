//
// Created by heuer on 24.05.17.
//

#include "Atom.h"

using namespace Eigen;
using namespace Elements;

Atom::Atom(Vector3d position, ElementType elementType)
        : Particle(position),
          elementType_(elementType)
{};

Atom::Atom(const Particle &particle, ElementType elementType)
        : Particle(particle),
          elementType_(elementType)
{}