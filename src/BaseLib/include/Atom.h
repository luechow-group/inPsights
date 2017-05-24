//
// Created by heuer on 23.05.17.
//

#ifndef AMOLQCGUI_ATOM_H
#define AMOLQCGUI_ATOM_H

#include "ElementTypes.h"
#include <Eigen/Core>

class Atom{
public:
    Atom(Elements::ElementType elementType, Eigen::Vector3d coordinates);

    Eigen::Vector3d coordinates() const { return coordinates_; };
    Elements::ElementType elementType() const { return elementType_; };

private:
    Elements::ElementType elementType_;
    Eigen::Vector3d coordinates_;
};


#endif //AMOLQCGUI_ATOM_H
