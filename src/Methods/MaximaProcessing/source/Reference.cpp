//
// Created by heuer on 11.04.19.
//


#include <Reference.h>
#include <Sample.h>
#include <Group.h>

Reference::Reference(double negLogSqrdProbabilityDensity, ElectronsVector maximum, size_t id)
        :
        negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity),
        maximum_(std::move(maximum)),
        sampleIds_({id})
{}

size_t Reference::ownId() const {
    return sampleIds_[0];
}

void Reference::mergeReference(Group::iterator &it) {
    assert((*it).representative()->ownId() != ownId() && "Self-associations are not allowed");

    sampleIds_.insert(
            sampleIds_.end(),
            make_move_iterator((*it).representative()->sampleIds_.begin()),
            make_move_iterator((*it).representative()->sampleIds_.end()));
}

void Reference::permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples){
    maximum_.permute(perm);

    for(const auto & i : sampleIds_){
        samples[i].permute(perm);
    }
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
