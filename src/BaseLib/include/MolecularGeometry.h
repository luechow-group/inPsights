//
// Created by Michael Heuer on 08.05.18.
//

#ifndef INPSIGHTS_MOLECULARGEOMETRY_H
#define INPSIGHTS_MOLECULARGEOMETRY_H

#include "ParticlesVector.h"
#include <list>

class MolecularGeometry{
public:
    MolecularGeometry();
    MolecularGeometry(AtomsVector atoms, ElectronsVector electrons);

    const AtomsVector& atoms() const;
    AtomsVector& atoms();

    const ElectronsVector& electrons() const;
    ElectronsVector & electrons();

    Particle<int> operator[](long i) const;

    std::pair<bool,long> findIndexByNumberedType(const NumberedType<int> &numberedType) const;

    NumberedType<int> findNumberedTypeByIndex(unsigned idx) const;

    long numberOfEntities() const;

    //TODO refactor Motifs::classifyMotifs since it is the only user of the method
    bool coreElectronQ(long i, double threshold = 0.01) const ;

    std::list<long> coreElectronsIndices(long k, double threshold = 0.01) const;

    std::list<long> coreElectronsIndices(double threshold = 0.01) const;

    std::list<long> valenceElectronsIndices(double threshold = 0.01) const;


    friend std::ostream& operator<<(std::ostream &os, const MolecularGeometry &mol) {
        os << mol.atoms() << std::endl;
        os << mol.electrons() << std::endl;
        return os;
    }

private:
    AtomsVector atoms_;
    ElectronsVector electrons_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<>
    struct convert<MolecularGeometry> {
        static Node encode(const MolecularGeometry& rhs);
        static bool decode(const Node& node, MolecularGeometry& rhs);
    };
    Emitter& operator<< (Emitter& out, const MolecularGeometry& p);
}


#endif //INPSIGHTS_MOLECULARGEOMETRY_H
