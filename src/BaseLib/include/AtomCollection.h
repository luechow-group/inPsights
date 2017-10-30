//
// Created by heuer on 24.05.17.
//

#ifndef AMOLQCGUI_ATOMCOLLECTION_H
#define AMOLQCGUI_ATOMCOLLECTION_H

#include <vector>
#include "ParticleCollection.h"
#include "ElementTypeCollection.h"
#include "Atom.h"


class AtomCollection : public ParticleCollection,public ElementTypeCollection{
public:
    AtomCollection() = default;
    explicit AtomCollection(const Eigen::VectorXd& positions);
    explicit AtomCollection(const VectorXd& positions, const VectorXi& spinTypes);

    Atom atom(long i);

    void insert (const Atom& atom, long i);
    void append (const Atom& atom);
    void prepend(const Atom& atom);


    void addAtom(const Elements::ElementType& elementType,
                 const double x, const double y, const double z);
    
    void addAtom(const Elements::ElementType& elementType,
                 const Vector3d& position );
};

#endif //AMOLQCGUI_ATOMCOLLECTION_H
