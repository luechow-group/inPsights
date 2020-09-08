// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Cluster.h>
#include <ParticleSelection.h>
#include <BestMatch.h>
#include <Reference.h>
#include <Eigen/Core>

std::string ToString::clusterToString(const Cluster &cluster) {
    std::stringstream ss;
    ss << cluster;
    return ss.str();
}

Cluster::Cluster()
        : std::vector<Cluster>(0),
          representative_(nullptr),
          selectedElectronsCount_(0){
}

Cluster::Cluster(Reference reference)
        : std::vector<Cluster>(0),
        representative_(std::make_shared<Reference>(std::move(reference))),
        selectedElectronsCount_(representative_->maximum().numberOfEntities()){
}

Cluster::Cluster(std::vector<Cluster>::size_type size)
        : std::vector<Cluster>(size),
        representative_(nullptr),
        selectedElectronsCount_(0){
}

Cluster::Cluster(std::initializer_list<Cluster> cluster)
        : std::vector<Cluster>(cluster),
        representative_( empty()? nullptr : front().representative()),
        selectedElectronsCount_(cluster.begin()->getSelectedElectronsCount()){
}

// Sort only this cluster and update the representative structure.
void Cluster::sort() {
    if(!isLeaf()) {
        std::sort(begin(), end());
        updateRepresentative();
    }
}

/* Sort all subclusters representative Reference objects by recursing down to the leaf level
 * and update all representative structures */
void Cluster::sortAll() {

    // sort all subclusters
    if(!isLeaf()) {
        for (auto &subcluster : *this)
            subcluster.sortAll();
    }

    // sort this cluster
    sort();
}

void Cluster::updateRepresentative() {
    representative_ = front().representative();
}

bool Cluster::isLeaf() const {
    return empty();
}

Cluster::size_type Cluster::numberOfLeaves() const {
    if(isLeaf()) {
        return 1;
    } else {
        Cluster::size_type numberOfLeaves = 0;
        for (auto &i : *this)
            numberOfLeaves += i.numberOfLeaves();
        return numberOfLeaves;
    }
}

void Cluster::permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic> &perm, std::vector<Sample>& samples) {
    if(isLeaf()) {
        representative()->permute(perm, samples);
    } else {
        for (auto &i : *this)
            i.permuteAll(perm, samples);
    }
}

void Cluster::permuteAll(const MolecularGeometry::Permutation &molecularPerm, std::vector<Sample> &samples) {
    if(isLeaf()) {
        representative()->permute(molecularPerm, samples);
    } else {
        for (auto &i : *this)
            i.permuteAll(molecularPerm, samples);
    }
}


ElectronsVector Cluster::electronsVectorFromAveragedPositionsVector(const AveragedPositionsVector & averagedPositionsVector) const {
    return {averagedPositionsVector.positions, representative()->maximum().typesVector()};
}

Cluster::AveragedPositionsVector Cluster::averagedSamplePositionsVector(const std::vector<Sample>& samples) const {
    auto sampleIds = allSampleIds();
    
    Eigen::VectorXd average = Eigen::VectorXd::Zero(
            representative()->maximum().numberOfEntities()
            *representative()->maximum().positionsVector().entityLength());
    for (const auto &sampleId : sampleIds) {
        average += samples[sampleId].sample_.positionsVector().asEigenVector();
    }
    average /= sampleIds.size();
    
    return {PositionsVector(average), static_cast<unsigned>(sampleIds.size())};
}

std::shared_ptr<Reference> Cluster::representative() {
    if (!isLeaf())
        return front().representative();
    else
        return representative_;
}

std::shared_ptr<const Reference> Cluster::representative() const {
    return std::const_pointer_cast<const Reference>(representative_);
}

bool Cluster::operator<(const Cluster &other) const {
    return *representative() < *other.representative();
}

Cluster &Cluster::operator+=(const Cluster &other) {

    if(isLeaf() && representative() != nullptr)
        emplace_back(*representative());

    if(other.isLeaf())
        emplace_back(other);
    else
        insert(end(), other.begin(), other.end());

    updateRepresentative();
    return *this;
}

void Cluster::makeSubcluster(std::vector<Cluster::iterator> its) {
    // sort the iterator list
    std::sort(its.begin(), its.end());

    Cluster subcluster;
    for (auto it : its)
        subcluster.emplace_back(*it);

    // reverse erase the iterators
    for(auto it = its.rbegin(); it != its.rend(); it++)
        erase(*it);

    emplace_back(subcluster);
    updateRepresentative();
}

std::vector<size_t> Cluster::allSampleIds() const {
    if(isLeaf()) {
        if(representative())
            return representative()->sampleIds();
        else
            return {};
    } else {
        std::vector<size_t> ids;
        for (const auto & subcluster : *this) {
            auto subclusterSampleIds = subcluster.allSampleIds();
            ids.insert(ids.end(), subclusterSampleIds.begin(), subclusterSampleIds.end());
        }
        return ids;
    }
}

std::ostream &operator<<(std::ostream &os, const Cluster &g) {
    os << "{";
    if(g.isLeaf()) {
        auto ids = g.allSampleIds();
        if(!ids.empty()) {
            for (auto it = ids.begin(); it != std::prev(ids.end()); it++)
                os << *it << ",";
            os << *std::prev(ids.end());
        }
    } else {
        for (auto it = g.begin(); it != std::prev(g.end()); it++)
            os << *it << ",";
        os << *std::prev(g.end());
    }
    os  << "}";
    
    return os;
}

long Cluster::getSelectedElectronsCount() const{
    return selectedElectronsCount_;
};

void Cluster::setSelectedElectronsCount(const long &count){
    selectedElectronsCount_ = count;
    for(auto & subcluster : (*this))
        subcluster.setSelectedElectronsCount(count);
};

void Cluster::permuteRelevantElectronsToFront(std::vector<Sample> & samples){
    Eigen::PermutationMatrix<Eigen::Dynamic> permutation;
    auto electronsNumber = (*this).representative()->maximum().numberOfEntities();

    for (auto & subCluster : *this) {
        auto subIndices = ParticleSelection::getRelevantIndices(subCluster.representative()->maximum());

        // permute all relevant electrons to the front
        subCluster.setSelectedElectronsCount(subIndices.size());
        permutation = BestMatch::getPermutationToFront(subIndices, electronsNumber);
        subCluster.permuteAll(permutation, samples);
    }
}