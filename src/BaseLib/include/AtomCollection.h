//
// Created by heuer on 24.05.17.
//

#ifndef AMOLQCGUI_ATOMCOLLECTION_H
#define AMOLQCGUI_ATOMCOLLECTION_H

#include "ParticleCollection.h"
#include "ElementTypeCollection.h"
#include "Atom.h"

class AtomCollection : public ParticleCollection{
public:
    AtomCollection() = default;
    explicit AtomCollection(const Eigen::VectorXd& positions);
    AtomCollection(const Eigen::VectorXd& positions, const Eigen::VectorXi& elementTypes);
    AtomCollection(const PositionsVector &positionsVector,
                   const ElementTypeCollection &elementTypeCollection);

    Atom operator[](long i) const;

    void insert (const Atom& atom, long i);
    void append (const Atom& atom);
    void prepend(const Atom& atom);
    void permute(long i, long j);


    const ElementTypeCollection& elementTypeCollection() const;
    ElementTypeCollection& elementTypeCollection();

    friend std::ostream& operator<<(std::ostream& os, const AtomCollection& ac);

private:
    ElementTypeCollection elementTypeCollection_;
};

#endif //AMOLQCGUI_ATOMCOLLECTION_H
