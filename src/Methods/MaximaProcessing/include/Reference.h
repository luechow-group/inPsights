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

#ifndef INPSIGHTS_REFERENCE_H
#define INPSIGHTS_REFERENCE_H

#include <ParticlesVector.h>
#include <Statistics.h>
#include <Group.h>
#include <MolecularSpectrum.h>

class Sample;

class Reference{
public:
    Reference(double negLogSqrdProbabilityDensity = std::numeric_limits<double>::max(),
              ElectronsVector maximum = {}, size_t id = 0);

    size_t ownId() const;

    //TODO is this the task of a container?
    void mergeReference(Group::iterator &it);

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples);

    bool operator<(const Reference& rhs) const;

    unsigned long count() const;

    double value() const;

    ElectronsVector maximum() const;

    std::vector<size_t> sampleIds() const;

    const SOAP::MolecularSpectrum& spectrum() const;
    void setSpectrum(SOAP::MolecularSpectrum spectrum);

private:
    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    std::vector<size_t> sampleIds_; //TODO use std::list or std::set?
    SOAP::MolecularSpectrum spectrum_;
};



#endif //INPSIGHTS_REFERENCE_H

