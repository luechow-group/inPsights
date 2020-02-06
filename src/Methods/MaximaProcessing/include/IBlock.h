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

#ifndef INPSIGHTS_IBLOCK_H
#define INPSIGHTS_IBLOCK_H

#include <Group.h>

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
    virtual void cluster(Group& group) = 0;
};

class IAnalyzer : public IBlock {
public:
    virtual void analyze(const Group& group) = 0;
};


#endif //INPSIGHTS_IBLOCK_H
