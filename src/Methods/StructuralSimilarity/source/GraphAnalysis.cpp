//
// Created by Michael Heuer on 2019-02-02.
//

#include <GraphAnalysis.h>

namespace GraphAnalysis {
    std::list<Eigen::Index> findConnectedVertices(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index startVertex) {
        assert(adjacencyMatrix.rows() == adjacencyMatrix.cols());

        Eigen::Index vertexCount = adjacencyMatrix.rows();
        assert(vertexCount > 0);
        assert(startVertex < vertexCount);

        Eigen::VectorXb marks = Eigen::VectorXb::Constant(vertexCount, 1, false);
        marks(startVertex) = true;


        std::queue<Eigen::Index> queue;
        queue.push(startVertex);

        // breadth first search
        while (!queue.empty()) {
            for (Eigen::Index i = 0; i < vertexCount; ++i)
                if (adjacencyMatrix(queue.front(), i) && !marks(i)) {
                    marks(i) = true;
                    queue.push(i);
                }
            queue.pop();
        }

        std::list<Eigen::Index> connectedVertices;

        for (Eigen::Index i = 0; i < vertexCount; ++i)
            if (marks(i)) connectedVertices.push_back(i);

        return connectedVertices;
    }

    std::vector<std::list<Eigen::Index>> findGraphClusters(const Eigen::MatrixXb &adjacencyMatrix) {
        assert(adjacencyMatrix.rows() == adjacencyMatrix.cols());

        Eigen::Index vertexCount = adjacencyMatrix.rows();
        assert(vertexCount > 0);

        std::list<Eigen::Index> remainingVertices;
        for (Eigen::Index i = 0; i < vertexCount; ++i)
            remainingVertices.push_back(i);

        std::vector<std::list<Eigen::Index>> clusters;

        do {
            auto connectedVertices = findConnectedVertices(adjacencyMatrix, remainingVertices.front());
            clusters.push_back(connectedVertices);

            std::list<Eigen::Index> difference;

            std::set_difference(
                    remainingVertices.begin(), remainingVertices.end(),
                    connectedVertices.begin(), connectedVertices.end(), std::back_inserter(difference));

            remainingVertices = difference;

        } while (!remainingVertices.empty());

        return clusters;
    }
}