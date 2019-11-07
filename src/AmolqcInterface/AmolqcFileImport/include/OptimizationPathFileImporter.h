/* Copyright (C) 2017-2019 Michael Heuer.
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

#ifndef INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H
#define INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H

#include "AmolqcImporter.h"
#include <ParticlesVectorCollection.h>

class PathElementDataEntry{
    PathElementDataEntry(unsigned long startLineIdx, unsigned long numberOfPathElements)
            : startLineIdx_(startLineIdx),
              numberOfPathElements_(numberOfPathElements){};
public:
    unsigned long startLineIdx_, numberOfPathElements_;
};

class OptimizationPathFileImporter : public AmolqcImporter{
public:
    explicit OptimizationPathFileImporter(const std::string& filename);

    ElectronsVectorCollection getPath(unsigned long k) const;

    AtomsVector getAtomsVector() const;

    PositionsVector getPositionsVector(unsigned long k, unsigned long m) const;

    unsigned long getNumberOfPaths() const;

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    AtomsVector atoms_;
    unsigned long numberOfNuclei_,beginOfElectronPositionBlocks_, numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfPaths_;
    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H
