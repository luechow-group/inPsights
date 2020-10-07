// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MAXIMAPROCESSOR_H
#define INPSIGHTS_MAXIMAPROCESSOR_H

#include "Maximum.h"
#include "Sample.h"
#include <Motifs.h>
#include <Statistics.h>
#include <VoxelCube.h>
#include <SOAPClusterer.h>
#include <EnergyResultsBundle.h>

namespace MotifEnergyCalculator {
    struct Result{
        VectorStatistics intraEnergies;
        TriangularMatrixStatistics interEnergies;
    };

    Result partition(const Cluster &cluster, const std::vector<Sample> &samples, const Motifs& motifs);

    void partitionLowerLevels(const Cluster &cluster,
                              const std::vector<Sample> &samples,
                              const Motifs& motifs,
                              VectorStatistics& intraMotifEnergyStats,
                              TriangularMatrixStatistics& interMotifEnergyStats);

    void partitionLowestLevel(const Cluster &cluster,
                              const std::vector<Sample> &samples,
                              const Motifs& motifs,
                              VectorStatistics& intraMotifEnergyStats,
                              TriangularMatrixStatistics& interMotifEnergyStats);

};

class MaximaProcessor {
public:

    MaximaProcessor(YAML::Emitter &yamlDocument, const std::vector<Sample> &samples, const AtomsVector &atoms);

    size_t addMaximum(const Maximum &maximum);

    size_t addAllMaxima(const Cluster &cluster);

    std::vector<ElectronsVector> getAllRepresentativeMaxima(const Cluster &cluster);

    void calculateStatistics(const Cluster &maxima,
                             const std::vector<std::vector<std::vector<size_t>>> &nucleiMergeLists,
                             const std::vector<size_t> &nucleiIndices,
                             const std::vector<DynamicMolecularSelection>& selections
                             );

    YAML::Node getYamlNode();

    std::string getYamlDocumentString();

private:
    YAML::Emitter &yamlDocument_;
    const std::vector<Sample> &samples_;
    AtomsVector atoms_;
    SingleValueStatistics valueStats_, EtotalStats_;
    VectorStatistics TeStats_, EeStats_, EnStats_;
    TriangularMatrixStatistics SeeStats_, VeeStats_, VnnStats_, ReeStats_;
    MatrixStatistics VenStats_, RenStats_;

    Eigen::MatrixXd Vnn_;
};

#endif //INPSIGHTS_MAXIMAPROCESSOR_H
