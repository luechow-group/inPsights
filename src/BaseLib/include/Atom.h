//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_ATOM_H
#define AMOLQCGUI_ATOM_H

#include "Particle.h"
#include "ElementType.h"

class Atom : public Particle{
public:
    Atom(const Eigen::Vector3d& position, const Elements::ElementType& elementType = Elements::ElementType::none);
    Atom(double x, double y, double z, const Elements::ElementType& elementType = Elements::ElementType::none);
    Atom(const Particle& particle, const Elements::ElementType& elementType = Elements::ElementType::none);

    Elements::ElementType elementType() const;
    void setElementType(const Elements::ElementType & elementType);

private:
    Elements::ElementType elementType_;
};


#endif //AMOLQCGUI_ATOM_H
