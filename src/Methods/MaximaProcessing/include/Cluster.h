// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTER_H
#define INPSIGHTS_CLUSTER_H

#include <vector>
#include <Eigen/Core>
#include <memory>
#include <Enumerate.h>
#include <Sample.h>
#include <MolecularGeometry.h>

// In order to use Cluster class,
// <Maximum.h> has to be included as well, due to forward declaration

class Cluster;
class Maximum;

namespace ToString {
    std::string clusterToString(const Cluster& cluster);
}

class Cluster : public std::vector<Cluster> {
public:
    Cluster();
    Cluster(const Cluster& group) = default;

    explicit Cluster(std::vector<Cluster>::size_type size);
    Cluster(std::initializer_list<Cluster> group);
    explicit Cluster(Maximum maximum);

    bool isLeaf() const;

    Cluster::size_type numberOfLeaves() const;

    void sort();
    void sortAll();
    void permuteRelevantElectronsToFront(std::vector<Sample> & samples);

    void getMaxima(std::vector<ElectronsVector> & maxima, std::size_t maximalNumber = std::numeric_limits<std::size_t>::max()) const;

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
    
    std::shared_ptr<const Maximum> representative() const;
    std::shared_ptr<Maximum> representative();

    bool operator<(const Cluster& other) const;

    std::vector<size_t > allSampleIds() const;

    friend std::ostream& operator<<(std::ostream& os, const Cluster & g);

    long getSelectedElectronsCount() const;
    void setSelectedElectronsCount(const long &count);

private:
    virtual void updateRepresentative();
    std::shared_ptr<Maximum> representative_;
    long selectedElectronsCount_;
};

inline Cluster operator+ (Cluster lhs, const Cluster& rhs) {
    lhs += rhs;
    return lhs;
}

#endif //INPSIGHTS_CLUSTER_H
