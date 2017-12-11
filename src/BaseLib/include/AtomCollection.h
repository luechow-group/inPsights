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
    AtomCollection(const Eigen::VectorXd& positions, const Eigen::VectorXi& elementTypes);
    AtomCollection(const ParticleCollection& particleCollection,
                   const ElementTypeCollection& elementTypeCollection);

    Atom atom(long i);

    void insert (const Atom& atom, long i);
    void append (const Atom& atom);
    void prepend(const Atom& atom);
    void permute(long i, long j);


    void addAtom(double x, double y, double z,
                 const Elements::ElementType &elementType = Elements::ElementType::none);
    
    void addAtom(const Eigen::Vector3d &position, const Elements::ElementType &elementType);
};

#endif //AMOLQCGUI_ATOMCOLLECTION_H
