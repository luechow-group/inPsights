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

#ifndef INPSIGHTS_CLUSTER_H
#define INPSIGHTS_CLUSTER_H

#include <vector>
#include <Eigen/Core>
#include <memory>
#include <Enumerate.h>
#include <Sample.h>
#include <MolecularGeometry.h>

// In order to use Cluster class,
// <Reference.h> has to be included as well, due to forward declaration

class Cluster;
class Reference;

namespace ToString {
    std::string clusterToString(const Cluster& cluster);
}

class Cluster : public std::vector<Cluster> {
public:
    Cluster();
    Cluster(const Cluster& group) = default;

    explicit Cluster(std::vector<Cluster>::size_type size);
    Cluster(std::initializer_list<Cluster> group);
    explicit Cluster(Reference reference);

    bool isLeaf() const;

    Cluster::size_type numberOfLeaves() const;

    void sort();
    void sortAll();
    void permuteRelevantElectronsToFront(std::vector<Sample> & samples);

    Cluster& operator+= (const Cluster& other);

    void makeSubcluster(std::vector<Cluster::iterator> its);

    void permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples);
    void permuteAll(const MolecularGeometry::Permutation& molecularPerm, std::vector<Sample>& samples);

    struct AveragedPositionsVector {
        PositionsVector positions;
        unsigned weight;
    };

    ElectronsVector electronsVectorFromAveragedPositionsVector(const AveragedPositionsVector& averagedPositionsVector)  const;

    AveragedPositionsVector averagedSamplePositionsVector(const std::vector<Sample>& samples) const;
    
    std::shared_ptr<const Reference> representative() const;
    std::shared_ptr<Reference> representative();

    bool operator<(const Cluster& other) const;

    std::vector<size_t > allSampleIds() const;

    friend std::ostream& operator<<(std::ostream& os, const Cluster & g);

    long getSelectedElectronsCount() const;
    void setSelectedElectronsCount(const long &count);

private:
    virtual void updateRepresentative();
    std::shared_ptr<Reference> representative_;
    long selectedElectronsCount_;
};

inline Cluster operator+ (Cluster lhs, const Cluster& rhs) {
    lhs += rhs;
    return lhs;
}

#endif //INPSIGHTS_CLUSTER_H
