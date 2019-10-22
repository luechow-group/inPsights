/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */
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
