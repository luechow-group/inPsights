// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_IPROCESS_H
#define INPSIGHTS_IPROCESS_H

#include <Cluster.h>

class IProcess{
public:
    enum class ProcessType {
        // Clusterer Processes
        IdentityClusterer,
        DistanceClusterer,
        DensityBasedClusterer,
        SOAPClusterer,
        ReferencePositionsClusterer,

        // Analyzer Processes
        ClusterNumberAnalyzer,
        TotalWeightDifferenceAnalyzer,

        // Unknown Process
        invalid
    };

    static ProcessType typeFromString(const std::string& processName);
};

class IClusterer : public IProcess {
public:
    explicit IClusterer(std::vector<Sample> & samples);

    virtual void cluster(Cluster& cluster) = 0;

protected:
    std::vector<Sample> &samples_;
};

class IAnalyzer : public IProcess {
public:
    virtual void analyze(const Cluster& cluster) = 0;
};


#endif //INPSIGHTS_IPROCESS_H
