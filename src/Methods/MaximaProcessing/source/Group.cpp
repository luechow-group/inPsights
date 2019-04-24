/* Copyright (C) 2019 Michael Heuer.
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

#include <Group.h>
#include <Reference.h>
#include <Eigen/Core>

Group::Group(Reference reference)
        : std::vector<Group>(0),
        representative_(std::make_shared<Reference>(std::move(reference))) {
}

Group::Group(std::vector<Group>::size_type size)
        : std::vector<Group>(size),
        representative_(nullptr){
}

Group::Group(std::initializer_list<Group> group)
        : std::vector<Group>(group),
        representative_( empty()? nullptr : front().representative()) {
}

// Sort only this group and update the representative structure.
void Group::sort() {
    if(!isLeaf()) {
        std::sort(begin(), end());
        updateRepresentative();
    }
}

/* Sort all subgroups representative Reference objects by recursing down to the leaf level
 * and update all representative structures */
void Group::sortAll() {

    // sort all subgroups
    if(!isLeaf()) {
        for (auto &subgroup : *this)
            subgroup.sortAll();
    }

    // sort this group
    sort();
}

void Group::updateRepresentative() {
    representative_ = front().representative();
}

bool Group::isLeaf() const {
    return empty();
}

Group::size_type Group::numberOfLeaves() const {
    if(isLeaf()) {
        return 1;
    } else {
        Group::size_type numberOfLeaves = 0;
        for (auto &i : *this)
            numberOfLeaves += i.numberOfLeaves();
        return numberOfLeaves;
    }
}

void Group::permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic> &perm, std::vector<Sample>& samples) {
    if(isLeaf()) {
        representative()->permute(perm, samples);
    } else {
        for (auto &i : *this)
            i.permuteAll(perm, samples);
    }
}

Group::AveragedPositionsVector Group::averagedPositionsVector() const {
    if (isLeaf())
        return {representative()->maximum().positionsVector(), 1};
    else {
        unsigned weight = 0;
        Eigen::VectorXd average = Eigen::VectorXd::Zero(
                representative()->maximum().numberOfEntities()
                *representative()->maximum().positionsVector().entityLength());
        for (const auto &subgroup : *this) {
            auto subgroupAverage = subgroup.averagedPositionsVector();
            average += double(subgroupAverage.weight)* subgroupAverage.positions.asEigenVector();
            weight += subgroupAverage.weight;
        }
        average /= weight;
        return {PositionsVector(average), weight};
    }
}

ElectronsVector Group::averagedRepresentativeElectronsVector() const {
    return {averagedPositionsVector().positions,representative()->maximum().typesVector()};
}

std::shared_ptr<Reference> Group::representative() {
    if (!isLeaf())
        return front().representative();
    else
        return representative_;
}

std::shared_ptr<const Reference> Group::representative() const {
    return std::const_pointer_cast<const Reference>(representative_);
}

bool Group::operator<(const Group &other) const {
    return *representative().get() < *other.representative();
}

Group &Group::operator+=(const Group &other) {

    if(isLeaf() && representative() != nullptr)
        emplace_back(*representative());

    if(other.isLeaf())
        emplace_back(other);
    else
        insert(end(), other.begin(), other.end());

    updateRepresentative();
    return *this;
}

void Group::makeSubgroup(std::vector<Group::iterator> its) {
    // sort the iterator list
    std::sort(its.begin(), its.end());

    Group subgroup;
    for (auto it : its)
        subgroup.emplace_back(*it);

    // reverse erase the iterators
    for(auto it = its.rbegin(); it != its.rend(); it++)
        erase(*it);

    emplace_back(subgroup);
    updateRepresentative();
}

std::vector<size_t> Group::allSampleIds() const {
    if(isLeaf())
        return representative()->sampleIds();
    else {
        std::vector<size_t> ids;
        for (auto subgroup : *this) {
            auto subgroupSampleIds = subgroup.allSampleIds();
            ids.insert(ids.end(), subgroupSampleIds.begin(), subgroupSampleIds.end());
        }
        return ids;
    }
}

std::ostream &operator<<(std::ostream &os, const Group &g) {
    os << "{";
    if(g.isLeaf()) {
        auto ids = g.allSampleIds();
        for (auto it = ids.begin(); it != std::prev(ids.end()); it++)
            os << *it << ",";
        os << *std::prev(ids.end());
    } else {
        for (auto it = g.begin(); it != std::prev(g.end()); it++)
            os << *it << ",";
        os << *std::prev(g.end());
    }
    os  << "}";
    
    return os;
}
