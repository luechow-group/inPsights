// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IBlock.h"

IBlock::BlockType IBlock::typeFromString(const std::string& blockName) {
    // Clusterer Blocks
    if(blockName == "IdentityClusterer")
        return BlockType::IdentityClusterer;
    else if(blockName == "PreClusterer")
        return BlockType::DistanceClusterer;
    else if(blockName == "DensityBasedClusterer")
        return BlockType::DensityBasedClusterer;
    else if(blockName == "ReferencePositionsClusterer")
        return BlockType::ReferencePositionsClusterer;
    else if(blockName == "SOAPClusterer")
        return BlockType::SOAPClusterer;

    // Analyzer Blocks
    else if(blockName == "ClusterNumberAnalyzer")
        return BlockType::ClusterNumberAnalyzer;
    else if(blockName == "TotalWeightDifferenceAnalyzer")
        return BlockType::TotalWeightDifferenceAnalyzer;

    // Unknown Block name
    else
        return BlockType::invalid;
};

IClusterer::IClusterer(std::vector<Sample> & samples)
: samples_(samples) {}
