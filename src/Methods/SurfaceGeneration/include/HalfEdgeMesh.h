// SPDX-License-Identifier: CC0-1.0
// from https://github.com/akuukka/quickhull

#include "Mesh.h"

#ifndef HALFEDGEMESH_H
#define HALFEDGEMESH_H

namespace quickhull {

    template<typename T, typename IndexType>
    class HalfEdgeMesh {
    public:

        struct HalfEdge {
            IndexType endVertex_;
            IndexType opp_;
            IndexType face_;
            IndexType next_;
        };

        struct Face {
            IndexType halfEdgeIndex_; // Index of one of the half edges of this face
        };

        std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> vertices_;
        std::vector<Face> faces_;
        std::vector<HalfEdge> halfEdges_;

        HalfEdgeMesh(const MeshBuilder <T> &builderObject, const VertexDataSource <T> &vertexData) {
            std::unordered_map<IndexType, IndexType> faceMapping;
            std::unordered_map<IndexType, IndexType> halfEdgeMapping;
            std::unordered_map<IndexType, IndexType> vertexMapping;

            size_t i = 0;
            for (const auto &face : builderObject.faces_) {
                if (!face.isDisabled()) {
                    faces_.push_back({static_cast<IndexType>(face.he_)});
                    faceMapping[i] = faces_.size() - 1;

                    const auto heIndices = builderObject.getHalfEdgeIndicesOfFace(face);
                    for (const auto heIndex : heIndices) {
                        const IndexType vertexIndex = builderObject.halfEdges_[heIndex].endVertex_;
                        if (vertexMapping.count(vertexIndex) == 0) {
                            vertices_.push_back(vertexData[vertexIndex]);
                            vertexMapping[vertexIndex] = vertices_.size() - 1;
                        }
                    }
                }
                i++;
            }

            i = 0;
            for (const auto &halfEdge : builderObject.halfEdges_) {
                if (!halfEdge.isDisabled()) {
                    halfEdges_.push_back({
                                                 static_cast<IndexType>(halfEdge.endVertex_),
                                                 static_cast<IndexType>(halfEdge.opp_),
                                                 static_cast<IndexType>(halfEdge.face_),
                                                 static_cast<IndexType>(halfEdge.next_)});
                    halfEdgeMapping[i] = halfEdges_.size() - 1;
                }
                i++;
            }

            for (auto &face : faces_) {
                assert(halfEdgeMapping.count(face.halfEdgeIndex_) == 1);
                face.halfEdgeIndex_ = halfEdgeMapping[face.halfEdgeIndex_];
            }

            for (auto &he : halfEdges_) {
                he.face_ = faceMapping[he.face_];
                he.opp_ = halfEdgeMapping[he.opp_];
                he.next_ = halfEdgeMapping[he.next_];
                he.endVertex_ = vertexMapping[he.endVertex_];
            }
        }

    };
}

#endif //HALFEDGEMESH_H
