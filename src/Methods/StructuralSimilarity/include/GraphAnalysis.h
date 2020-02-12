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

#ifndef INPSIGHTS_GRAPHANALYSIS_H
#define INPSIGHTS_GRAPHANALYSIS_H

#include <Eigen/Core>
#include <queue>
#include <list>
#include <algorithm>
#include <map>

namespace Eigen {
    using MatrixXb = Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXb = Matrix<bool, Eigen::Dynamic, 1>;
}
namespace GraphAnalysis {
    // filter to obtain an adjacency matrix
    Eigen::MatrixXb filter(const Eigen::MatrixXd & matrix, double threshold = 1.0);

    Eigen::MatrixXb lowerOrEqualFilter(const Eigen::MatrixXd & matrix, double threshold = 1.0);

    std::list<Eigen::Index> findConnectedVertices(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index startVertex);

    std::vector<std::list<Eigen::Index>> findGraphClusters(const Eigen::MatrixXb &adjacencyMatrix);

    std::map<std::size_t, std::size_t> findMergeMap(
            std::vector<std::list<Eigen::Index>> subsets,
            std::vector<std::list<Eigen::Index>> referenceSets);


    template<typename K, typename V>
    bool findByValue(std::vector<K> & vec, std::map<K, V> mapOfElement, V value)
    {
        bool bResult = false;
        auto it = mapOfElement.begin();
        // Iterate through the map
        while(it != mapOfElement.end())
        {
            // Check if value of this entry matches with given value
            if(it->second == value)
            {
                // Yes found
                bResult = true;
                // Push the key in given map
                vec.push_back(it->first);
            }
            // Go to next entry in map
            it++;
        }
        return bResult;
    }

}

#endif //INPSIGHTS_GRAPHANALYSIS_H
