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

#ifndef INPSIGHTS_MAXIMAPROCESSOR_H
#define INPSIGHTS_MAXIMAPROCESSOR_H

#include "Reference.h"
#include "Sample.h"
#include <Motifs.h>
#include <Statistics.h>
#include <VoxelCube.h>
#include <SOAPClusterer.h>

namespace MotifEnergyCalculator {
    struct Result{
        VectorStatistics intraEnergies;
        TriangularMatrixStatistics interEnergies;
    };

    Result partition(const Group &group, const std::vector<Sample> &samples, const Motifs& motifs);

    void partitionLowerLevels(const Group &group,
                              const std::vector<Sample> &samples,
                              const Motifs& motifs,
                              VectorStatistics& intraMotifEnergyStats,
                              TriangularMatrixStatistics& interMotifEnergyStats);

    void partitionLowestLevel(const Group &group,
                              const std::vector<Sample> &samples,
                              const Motifs& motifs,
                              VectorStatistics& intraMotifEnergyStats,
                              TriangularMatrixStatistics& interMotifEnergyStats);

};

class MaximaProcessor {
public:

    MaximaProcessor(YAML::Emitter &yamlDocument, const std::vector<Sample> &samples, const AtomsVector &atoms);

    size_t addReference(const Reference &reference);

    size_t addAllReferences(const Group &group);

    std::vector<ElectronsVector> getAllRepresentativeMaxima(const Group &group);

    void calculateStatistics(const Group &maxima,
                             const std::vector<std::vector<std::vector<size_t>>> &nucleiMergeLists);

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
