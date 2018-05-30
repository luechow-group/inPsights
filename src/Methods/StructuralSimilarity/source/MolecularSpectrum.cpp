//
// Created by Michael Heuer on 27.05.18.
//
#include "MolecularSpectrum.h"

MolecularSpectrum::MolecularSpectrum(MolecularGeometry molecule)
        : molecule_(molecule)
{
    assert(ParticleKit::isSubsetQ(molecule));
    NeighborhoodExpander expander;
    molecularCenters_ = expander.computeExpansions(molecule);
}
