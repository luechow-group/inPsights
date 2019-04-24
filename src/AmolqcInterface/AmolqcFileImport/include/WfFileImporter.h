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

#ifndef WFFILEIMPORTER_H
#define WFFILEIMPORTER_H

#include "AmolqcImporter.h"
#include <ParticlesVector.h>

enum JastrowTypes{
    none,
    IC,
    SM
};

class WfFileImporter : public AmolqcImporter {
public:
    explicit WfFileImporter(const std::string &filename);

    std::pair<bool, unsigned long> findTag(const std::string &tag,
                                           unsigned long startLine = 0);

    AtomsVector getAtomsVector();
    unsigned long getNumberOfElectrons();
    unsigned long getNumberOfAlphaElectrons();
    unsigned long getNumberOfBetaElectrons();
    long getCharge(){ return charge_; };
    unsigned long getMultiplicity(){ return multiplicity_; };
    std::string getBasisSet(){ return basis_; };

private:
    void readGeneralBlock();

    bool bohrQ_;
    long charge_;
    unsigned long numberOfNuclei_,multiplicity_;
    std::string basis_, title_, jastrow_;
};

#endif //WFFILEIMPORTER_H
