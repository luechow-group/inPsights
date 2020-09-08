// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later


#include <Reference.h>
#include <Sample.h>
#include <Cluster.h>

Reference::Reference(const AtomsVector &nuclei, double negLogSqrdProbabilityDensity, ElectronsVector maximum, size_t id)
        :
        negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity),
        maximum_(std::move(maximum)),
        sampleIds_({id}),
        nucleiRef(nuclei),
        nuclearPermutation_(nuclei.numberOfEntities())
{
    nuclearPermutation_.setIdentity();
}

size_t Reference::ownId() const {
    return sampleIds_[0];
}

void Reference::mergeReference(Cluster::iterator &it) {
    assert(it->representative()->ownId() != ownId() && "Self-associations are not allowed.");
    assert(it->representative()->nuclei() == nuclei() && "The nuclear geometries must match.");
    sampleIds_.insert(
            sampleIds_.end(),
            make_move_iterator((*it).representative()->sampleIds_.begin()),
            make_move_iterator((*it).representative()->sampleIds_.end()));
}

void Reference::permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples){
    maximum_.permute(perm);

    for(const auto & i : sampleIds_)
        samples[i].permute(perm);
}

void Reference::permute(const MolecularGeometry::Permutation& molPerm, std::vector<Sample>& samples){
    maximum_.permute(molPerm.electronicPermutation);

    for(const auto & i : sampleIds_)
        samples[i].permute(molPerm.electronicPermutation);

    nuclearPermutation_ = molPerm.nuclearPermutation * nuclearPermutation_;
}

bool Reference::operator<(const Reference& rhs) const {
    return negLogSqrdProbabilityDensity_< rhs.negLogSqrdProbabilityDensity_;
}

unsigned long Reference::count() const {
    return sampleIds_.size();
}

double Reference::value() const {
    return negLogSqrdProbabilityDensity_;
}

const AtomsVector &Reference::nuclei() const {
    return nucleiRef;
}

ElectronsVector Reference::maximum() const {
    return maximum_;
}

std::vector<size_t> Reference::sampleIds() const {
    return sampleIds_;
}

void Reference::setSpectrum(SOAP::MolecularSpectrum spectrum) {
    assert(spectrum.molecule_.electrons().typesVector() == maximum().typesVector()
    && "The TypesVectors must be identical.");
    spectrum_ = std::move(spectrum);
}

const SOAP::MolecularSpectrum& Reference::spectrum() const {
    return spectrum_;
}

Eigen::PermutationMatrix<Eigen::Dynamic> Reference::nuclearPermutation() const {
    assert( nuclearPermutation_.size() == nuclei().numberOfEntities());
    return nuclearPermutation_;
}
