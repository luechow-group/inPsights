//
// Created by Michael Heuer on 2019-02-02.
//

#ifndef INPSIGHTS_GRAPHANALYSIS_H
#define INPSIGHTS_GRAPHANALYSIS_H

#include <Eigen/Core>
#include <queue>
#include <list>
#include <algorithm>

namespace Eigen {
    using MatrixXb = Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXb = Matrix<bool, Eigen::Dynamic, 1>;
}
namespace GraphAnalysis {
    std::list<Eigen::Index> findConnectedVertices(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index startVertex);
    
    std::vector<std::list<Eigen::Index>> findGraphClusters(const Eigen::MatrixXb &adjacencyMatrix);
}


#endif //INPSIGHTS_GRAPHANALYSIS_H
