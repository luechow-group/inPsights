//
// Created by Michael Heuer on 02.10.18.
//

#ifndef INPSIGHTS_MAXIMAPROCESSOR_H
#define INPSIGHTS_MAXIMAPROCESSOR_H

#include "Reference.h"
#include "Sample.h"
#include <Motifs.h>
#include <Statistics.h>
#include <VoxelCube.h>
#include <SOAPClusterer.h>

class MaximaProcessor {
public:

    MaximaProcessor(YAML::Emitter& yamlDocument, const std::vector<Sample> &samples, AtomsVector atoms);

    size_t addReference(const Reference &reference);

    size_t addAllReferences(const Group &group);

    std::vector<ElectronsVector> getAllRepresentativeMaxima(const Group &group);

    void doMotifBasedEnergyPartitioning(const Group &group);

    void calculateStatistics(const Group &maxima);

    YAML::Node getYamlNode();

    std::string getYamlDocumentString();

private:
    YAML::Emitter& yamlDocument_;
    const std::vector<Sample> &samples_;
    AtomsVector atoms_;
    Motifs motifs_;
    SingleValueStatistics valueStats_, EtotalStats_;
    VectorStatistics TeStats_, EeStats_, EnStats_, intraMotifEnergyStats_;
    TriangularMatrixStatistics SeeStats_, VeeStats_, VnnStats_, interMotifEnergyStats_, ReeStats_;
    MatrixStatistics VenStats_, RenStats_;

    Eigen::MatrixXd Vnn_;
};

#endif //INPSIGHTS_MAXIMAPROCESSOR_H
