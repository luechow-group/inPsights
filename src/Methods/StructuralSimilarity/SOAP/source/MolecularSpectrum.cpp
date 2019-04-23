//
// Created by Michael Heuer on 27.05.18.
//
#include "MolecularSpectrum.h"

using namespace SOAP;

MolecularSpectrum::MolecularSpectrum(MolecularGeometry molecule)
        : molecule_(molecule)
{
    assert(molecule.numberOfEntities() > 0 && "The molecule must contain at least one particles.");
    assert(ParticleKit::isSubsetQ(molecule));
    NeighborhoodExpander expander;
    molecularCenters_ = expander.computeMolecularExpansions(molecule);
}
