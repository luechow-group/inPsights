// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later
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
