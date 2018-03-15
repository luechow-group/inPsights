//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCPP_ATOM_H
#define AMOLQCPP_ATOM_H

#include "Particle.h"
#include "ElementType.h"

class Atom : public Particle{
public:
    Atom(const Eigen::Vector3d& position, const Elements::ElementType& elementType = Elements::ElementType::none);
    Atom(double x, double y, double z, const Elements::ElementType& elementType = Elements::ElementType::none);
    Atom(const Particle& particle, const Elements::ElementType& elementType = Elements::ElementType::none);

    Elements::ElementType elementType() const;
    void setElementType(const Elements::ElementType & elementType);

    friend std::ostream& operator<< (std::ostream& os, const Atom& atom);

    std::string toString() const override;
    int charge() const override;

private:
    Elements::ElementType elementType_;
};


#endif //AMOLQCPP_ATOM_H
