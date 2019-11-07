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

#ifndef INPSIGHTS_REFFILEIMPORTER_H
#define INPSIGHTS_REFFILEIMPORTER_H

#include "AmolqcImporter.h"
#include <ParticlesVectorCollection.h>

class RefFileImporter : public AmolqcImporter{
public:
    RefFileImporter(const std::string& filename);

    AtomsVector getAtomsVector();

    SpinTypesVector getSpinTypesVector() const;
    PositionsVector getPositionsVector(unsigned long k, unsigned long m) const;
    ElectronsVector getMaximaStructure(unsigned long k, unsigned long m) const;
    ElectronsVectorCollection getAllSubstructures(unsigned long k) const;

    unsigned long getNumberOfMaxima(unsigned long k, unsigned long m) const;
    double getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const;

    unsigned long numberOfSuperstructures();

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    unsigned long numberOfNuclei_,
            numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfSuperstructures_,
            totalNumberOfMaxima_,
            maximalNumberOfSubstructures; ;
    const unsigned numberOfLinesAboveCoordinatesBlock = 2;

    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //INPSIGHTS_REFFILEIMPORTER_H
