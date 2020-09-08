// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_IBLOCK_H
#define INPSIGHTS_IBLOCK_H

#include <Cluster.h>

class IBlock{
public:
    enum class BlockType {
        // Clusterer Blocks
        IdentityClusterer,
        DistanceClusterer,
        DensityBasedClusterer,
        SOAPClusterer,
        ReferencePositionsClusterer,

        // Analyzer Blocks
        ClusterNumberAnalyzer,
        TotalWeightDifferenceAnalyzer,

        // Unknown Block name
        invalid
    };

    static BlockType typeFromString(const std::string& blockName);
};

class IClusterer : public IBlock {
public:
    explicit IClusterer(std::vector<Sample> & samples);

    virtual void cluster(Cluster& cluster) = 0;

protected:
    std::vector<Sample> &samples_;
};

class IAnalyzer : public IBlock {
public:
    virtual void analyze(const Cluster& cluster) = 0;
};


#endif //INPSIGHTS_IBLOCK_H
