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

#include "IBlock.h"

IBlock::BlockType IBlock::typeFromString(const std::string& blockName) {
    // Clusterer Blocks
    if(blockName == "IdentityClusterer")
        return BlockType::IdentityClusterer;
    else if(blockName == "DistanceClusterer")
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
