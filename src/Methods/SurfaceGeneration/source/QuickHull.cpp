#include "QuickHull.h"
#include "MathUtils.h"
#include "Types.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <deque>
#include <limits>
#include "Mesh.h"

namespace quickhull {

    template<>
    float defaultEps() {
        return 0.0001f;
    }

    template<>
    double defaultEps() {
        return 0.0000001;
    }

    /*
     * Implementation of the algorithm
     */

    template<typename T>
    ConvexHull<T> QuickHull<T>::getConvexHull(const std::vector<Eigen::Matrix<T, 3, 1>> &pointCloud, bool CCW,
                                              bool useOriginalIndices, T epsilon) {
        VertexDataSource<T> vertexDataSource(pointCloud);
        return getConvexHull(vertexDataSource, CCW, useOriginalIndices, epsilon);
    }

    template<typename T>
    ConvexHull<T> QuickHull<T>::getConvexHull(const Eigen::Matrix<T, 3, 1> *vertexData, size_t vertexCount, bool CCW,
                                              bool useOriginalIndices, T epsilon) {
        VertexDataSource<T> vertexDataSource(vertexData, vertexCount);
        return getConvexHull(vertexDataSource, CCW, useOriginalIndices, epsilon);
    }

    template<typename T>
    ConvexHull<T>
    QuickHull<T>::getConvexHull(const T *vertexData, size_t vertexCount, bool CCW, bool useOriginalIndices, T epsilon) {
        VertexDataSource<T> vertexDataSource((const Eigen::Matrix<T, 3, 1> *) vertexData, vertexCount);
        return getConvexHull(vertexDataSource, CCW, useOriginalIndices, epsilon);
    }

    template<typename T>
    HalfEdgeMesh<T, IndexType>
    QuickHull<T>::getConvexHullAsMesh(const T *vertexData, size_t vertexCount, bool CCW, T epsilon) {
        VertexDataSource<T> vertexDataSource((const Eigen::Matrix<T, 3, 1> *) vertexData, vertexCount);
        buildMesh(vertexDataSource, CCW, false, epsilon);
        return HalfEdgeMesh<T, IndexType>(mesh_, vertexData_);
    }

    template<typename T>
    void QuickHull<T>::buildMesh(const VertexDataSource<T> &pointCloud, bool CCW, bool useOriginalIndices, T epsilon) {
        if (pointCloud.size() == 0) {
            mesh_ = MeshBuilder<T>();
            return;
        }
        vertexData_ = pointCloud;

        // Very first: find extreme values and use them to compute the scale of the point cloud.
        extremeValues_ = getExtremeValues();
        scale_ = getScale(extremeValues_);

        // Epsilon we use depends on the scale
        epsilon_ = epsilon * scale_;
        epsilonSquared_ = epsilon_ * epsilon_;

        // Reset diagnostics
        diagnostics_ = DiagnosticsData();

        planar_ = false; // The planar case happens when all the points appear to lie on a two dimensional subspace of R^3.
        createConvexHalfEdgeMesh();
        if (planar_) {
            const size_t extraPointIndex = planarPointCloudTemp_.size() - 1;
            for (auto &he : mesh_.halfEdges_) {
                if (he.endVertex_ == extraPointIndex) {
                    he.endVertex_ = 0;
                }
            }
            vertexData_ = pointCloud;
            planarPointCloudTemp_.clear();
        }
    }

    template<typename T>
    ConvexHull<T>
    QuickHull<T>::getConvexHull(const VertexDataSource<T> &pointCloud, bool CCW, bool useOriginalIndices, T epsilon) {
        buildMesh(pointCloud, CCW, useOriginalIndices, epsilon);
        return ConvexHull<T>(mesh_, vertexData_, CCW, useOriginalIndices);
    }

    template<typename T>
    void QuickHull<T>::createConvexHalfEdgeMesh() {
        // Temporary variables used during iteration
        std::vector<IndexType> visibleFaces;
        std::vector<IndexType> horizonEdges;
        struct FaceData {
            IndexType faceIndex_;
            IndexType enteredFromHalfEdge_; // If the face turns out not to be visible, this half edge will be marked as horizon edge
            FaceData(IndexType fi, IndexType he) : faceIndex_(fi), enteredFromHalfEdge_(he) {}
        };
        std::vector<FaceData> possiblyVisibleFaces;

        // Compute base tetrahedron
        mesh_ = getInitialTetrahedron();
        assert(mesh_.faces_.size() == 4);

        // Init face stack with those faces that have points assigned to them
        std::deque<IndexType> faceList;
        for (size_t i = 0; i < 4; i++) {
            auto &f = mesh_.faces_[i];
            if (f.pointsOnPositiveSide_ && f.pointsOnPositiveSide_->size() > 0) {
                faceList.push_back(i);
                f.inFaceStack_ = 1;
            }
        }

        // Process faces until the face list is empty.
        size_t iter = 0;
        while (!faceList.empty()) {
            iter++;
            if (iter == std::numeric_limits<size_t>::max()) {
                // Visible face traversal marks visited faces with iteration counter (to mark that the face has been visited on this iteration) and the max value represents unvisited faces. At this point we have to reset iteration counter. This shouldn't be an
                // issue on 64 bit machines.
                iter = 0;
            }

            const IndexType topFaceIndex = faceList.front();
            faceList.pop_front();

            auto &tf = mesh_.faces_[topFaceIndex];
            tf.inFaceStack_ = 0;

            assert(!tf.pointsOnPositiveSide_ || tf.pointsOnPositiveSide_->size() > 0);
            if (!tf.pointsOnPositiveSide_ || tf.isDisabled()) {
                continue;
            }

            // Pick the most distant point to this triangle plane as the point to which we extrude
            const Eigen::Matrix<T, 3, 1> &activePoint = vertexData_[tf.mostDistantPoint_];
            const size_t activePointIndex = tf.mostDistantPoint_;

            // Find out the faces that have our active point on their positive side (these are the "visible faces"). The face on top of the stack of course is one of them. At the same time, we create a list of horizon edges.
            horizonEdges.clear();
            possiblyVisibleFaces.clear();
            visibleFaces.clear();
            possiblyVisibleFaces.emplace_back(topFaceIndex, std::numeric_limits<size_t>::max());
            while (possiblyVisibleFaces.size()) {
                const auto faceData = possiblyVisibleFaces.back();
                possiblyVisibleFaces.pop_back();
                auto &pvf = mesh_.faces_[faceData.faceIndex_];
                assert(!pvf.isDisabled());

                if (pvf.visibilityCheckedOnIteration_ == iter) {
                    if (pvf.isVisibleFaceOnCurrentIteration_) {
                        continue;
                    }
                } else {
                    const Plane<T> &P = pvf.P_;
                    pvf.visibilityCheckedOnIteration_ = iter;
                    const T d = P.N_.dot(activePoint) + P.D_;
                    if (d > 0) {
                        pvf.isVisibleFaceOnCurrentIteration_ = 1;
                        pvf.horizonEdgesOnCurrentIteration_ = 0;
                        visibleFaces.push_back(faceData.faceIndex_);
                        for (auto heIndex : mesh_.getHalfEdgeIndicesOfFace(pvf)) {
                            if (mesh_.halfEdges_[heIndex].opp_ != faceData.enteredFromHalfEdge_) {
                                possiblyVisibleFaces.emplace_back(
                                        mesh_.halfEdges_[mesh_.halfEdges_[heIndex].opp_].face_, heIndex);
                            }
                        }
                        continue;
                    }
                    assert(faceData.faceIndex_ != topFaceIndex);
                }

                // The face is not visible. Therefore, the halfedge we came from is part of the horizon edge.
                pvf.isVisibleFaceOnCurrentIteration_ = 0;
                horizonEdges.push_back(faceData.enteredFromHalfEdge_);
                // Store which half edge is the horizon edge. The other half edges of the face will not be part of the final mesh so their data slots can by recycled.
                const auto halfEdges = mesh_.getHalfEdgeIndicesOfFace(
                        mesh_.faces_[mesh_.halfEdges_[faceData.enteredFromHalfEdge_].face_]);
                const std::int8_t ind = (halfEdges[0] == faceData.enteredFromHalfEdge_) ? 0 : (halfEdges[1] ==
                                                                                               faceData.enteredFromHalfEdge_
                                                                                               ? 1 : 2);
                mesh_.faces_[mesh_.halfEdges_[faceData.enteredFromHalfEdge_].face_].horizonEdgesOnCurrentIteration_ |= (
                        1 << ind);
            }
            const size_t horizonEdgeCount = horizonEdges.size();

            // Order horizon edges so that they form a loop. This may fail due to numerical instability in which case we give up trying to solve horizon edge for this point and accept a minor degeneration in the convex hull.
            if (!reorderHorizonEdges(horizonEdges)) {
                diagnostics_.failedHorizonEdges_++;
                std::cerr << "Failed to solve horizon edge." << std::endl;
                auto it = std::find(tf.pointsOnPositiveSide_->begin(), tf.pointsOnPositiveSide_->end(),
                                    activePointIndex);
                tf.pointsOnPositiveSide_->erase(it);
                if (tf.pointsOnPositiveSide_->size() == 0) {
                    reclaimToIndexVectorPool(tf.pointsOnPositiveSide_);
                }
                continue;
            }

            // Except for the horizon edges, all half edges of the visible faces can be marked as disabled. Their data slots will be reused.
            // The faces will be disabled as well, but we need to remember the points that were on the positive side of them - therefore
            // we save pointers to them.
            newFaceIndices_.clear();
            newHalfEdgeIndices_.clear();
            disabledFacePointVectors_.clear();
            size_t disableCounter = 0;
            for (auto faceIndex : visibleFaces) {
                auto &disabledFace = mesh_.faces_[faceIndex];
                auto halfEdges = mesh_.getHalfEdgeIndicesOfFace(disabledFace);
                for (size_t j = 0; j < 3; j++) {
                    if ((disabledFace.horizonEdgesOnCurrentIteration_ & (1 << j)) == 0) {
                        if (disableCounter < horizonEdgeCount * 2) {
                            // Use on this iteration
                            newHalfEdgeIndices_.push_back(halfEdges[j]);
                            disableCounter++;
                        } else {
                            // Mark for reusal on later iteration step
                            mesh_.disableHalfEdge(halfEdges[j]);
                        }
                    }
                }
                // Disable the face, but retain pointer to the points that were on the positive side of it. We need to assign those points
                // to the new faces we create shortly.
                auto t = std::move(mesh_.disableFace(faceIndex));
                if (t) {
                    assert(t->size()); // Because we should not assign point vectors to faces unless needed...
                    disabledFacePointVectors_.push_back(std::move(t));
                }
            }
            if (disableCounter < horizonEdgeCount * 2) {
                const size_t newHalfEdgesNeeded = horizonEdgeCount * 2 - disableCounter;
                for (size_t i = 0; i < newHalfEdgesNeeded; i++) {
                    newHalfEdgeIndices_.push_back(mesh_.addHalfEdge());
                }
            }

            // Create new faces using the edgeloop
            for (size_t i = 0; i < horizonEdgeCount; i++) {
                const IndexType AB = horizonEdges[i];

                auto horizonEdgeVertexIndices = mesh_.getVertexIndicesOfHalfEdge(mesh_.halfEdges_[AB]);
                IndexType A, B, C;
                A = horizonEdgeVertexIndices[0];
                B = horizonEdgeVertexIndices[1];
                C = activePointIndex;

                const IndexType newFaceIndex = mesh_.addFace();
                newFaceIndices_.push_back(newFaceIndex);

                const IndexType CA = newHalfEdgeIndices_[2 * i + 0];
                const IndexType BC = newHalfEdgeIndices_[2 * i + 1];

                mesh_.halfEdges_[AB].next_ = BC;
                mesh_.halfEdges_[BC].next_ = CA;
                mesh_.halfEdges_[CA].next_ = AB;

                mesh_.halfEdges_[BC].face_ = newFaceIndex;
                mesh_.halfEdges_[CA].face_ = newFaceIndex;
                mesh_.halfEdges_[AB].face_ = newFaceIndex;

                mesh_.halfEdges_[CA].endVertex_ = A;
                mesh_.halfEdges_[BC].endVertex_ = C;

                auto &newFace = mesh_.faces_[newFaceIndex];

                const Eigen::Matrix<T, 3, 1> planeNormal = mathutils::getTriangleNormal(vertexData_[A], vertexData_[B],
                                                                                        activePoint);
                newFace.P_ = Plane<T>(planeNormal, activePoint);
                newFace.he_ = AB;

                mesh_.halfEdges_[CA].opp_ = newHalfEdgeIndices_[i > 0 ? i * 2 - 1 : 2 * horizonEdgeCount - 1];
                mesh_.halfEdges_[BC].opp_ = newHalfEdgeIndices_[((i + 1) * 2) % (horizonEdgeCount * 2)];
            }

            // Assign points that were on the positive side of the disabled faces to the new faces.
            for (auto &disabledPoints : disabledFacePointVectors_) {
                assert(disabledPoints);
                for (const auto &point : *(disabledPoints)) {
                    if (point == activePointIndex) {
                        continue;
                    }
                    for (size_t j = 0; j < horizonEdgeCount; j++) {
                        if (addPointToFace(mesh_.faces_[newFaceIndices_[j]], point)) {
                            break;
                        }
                    }
                }
                // The points are no longer needed: we can move them to the vector pool for reuse.
                reclaimToIndexVectorPool(disabledPoints);
            }

            // Increase face stack size if needed
            for (const auto newFaceIndex : newFaceIndices_) {
                auto &newFace = mesh_.faces_[newFaceIndex];
                if (newFace.pointsOnPositiveSide_) {
                    assert(newFace.pointsOnPositiveSide_->size() > 0);
                    if (!newFace.inFaceStack_) {
                        faceList.push_back(newFaceIndex);
                        newFace.inFaceStack_ = 1;
                    }
                }
            }
        }

        // Cleanup
        indexVectorPool_.clear();
    }

    /*
     * Private helper functions
     */

    template<typename T>
    std::array<IndexType, 6> QuickHull<T>::getExtremeValues() {
        std::array<IndexType, 6> outIndices{0, 0, 0, 0, 0, 0};
        T extremeVals[6] = {
                vertexData_[0].x(), vertexData_[0].x(), vertexData_[0].y(), vertexData_[0].y(), vertexData_[0].z(),
                vertexData_[0].z()};
        const size_t vCount = vertexData_.size();
        for (size_t i = 1; i < vCount; i++) {
            const Eigen::Matrix<T, 3, 1> &pos = vertexData_[i];
            if (pos.x() > extremeVals[0]) {
                extremeVals[0] = pos.x();
                outIndices[0] = (IndexType) i;
            } else if (pos.x() < extremeVals[1]) {
                extremeVals[1] = pos.x();
                outIndices[1] = (IndexType) i;
            }
            if (pos.y() > extremeVals[2]) {
                extremeVals[2] = pos.y();
                outIndices[2] = (IndexType) i;
            } else if (pos.y() < extremeVals[3]) {
                extremeVals[3] = pos.y();
                outIndices[3] = (IndexType) i;
            }
            if (pos.z() > extremeVals[4]) {
                extremeVals[4] = pos.z();
                outIndices[4] = (IndexType) i;
            } else if (pos.z() < extremeVals[5]) {
                extremeVals[5] = pos.z();
                outIndices[5] = (IndexType) i;
            }
        }
        return outIndices;
    }

    template<typename T>
    bool QuickHull<T>::reorderHorizonEdges(std::vector<IndexType> &horizonEdges) {
        const size_t horizonEdgeCount = horizonEdges.size();
        for (size_t i = 0; i < horizonEdgeCount - 1; i++) {
            const IndexType endVertex = mesh_.halfEdges_[horizonEdges[i]].endVertex_;
            bool foundNext = false;
            for (size_t j = i + 1; j < horizonEdgeCount; j++) {
                const IndexType beginVertex = mesh_.halfEdges_[mesh_.halfEdges_[horizonEdges[j]].opp_].endVertex_;
                if (beginVertex == endVertex) {
                    std::swap(horizonEdges[i + 1], horizonEdges[j]);
                    foundNext = true;
                    break;
                }
            }
            if (!foundNext) {
                return false;
            }
        }
        assert(mesh_.halfEdges_[horizonEdges[horizonEdges.size() - 1]].endVertex_ ==
               mesh_.halfEdges_[mesh_.halfEdges_[horizonEdges[0]].opp_].endVertex_);
        return true;
    }

    template<typename T>
    T QuickHull<T>::getScale(const std::array<IndexType, 6> &extremeValues) {
        T s = 0;
        for (size_t i = 0; i < 6; i++) {
            const T *v = (const T *) (&vertexData_[extremeValues[i]]);
            v += i / 2;
            auto a = std::abs(*v);
            if (a > s) {
                s = a;
            }
        }
        return s;
    }

    template<typename T>
    MeshBuilder<T> QuickHull<T>::getInitialTetrahedron() {
        const size_t vertexCount = vertexData_.size();

        // If we have at most 4 points, just return a degenerate tetrahedron:
        if (vertexCount <= 4) {
            IndexType v[4] = {0, std::min((size_t) 1, vertexCount - 1), std::min((size_t) 2, vertexCount - 1),
                              std::min((size_t) 3, vertexCount - 1)};
            const Eigen::Matrix<T, 3, 1> N = mathutils::getTriangleNormal(vertexData_[v[0]], vertexData_[v[1]],
                                                                          vertexData_[v[2]]);
            const Plane<T> trianglePlane(N, vertexData_[v[0]]);
            if (trianglePlane.isPointOnPositiveSide(vertexData_[v[3]])) {
                std::swap(v[0], v[1]);
            }
            return MeshBuilder<T>(v[0], v[1], v[2], v[3]);
        }

        // Find two most distant extreme points.
        T maxD = epsilonSquared_;
        std::pair<IndexType, IndexType> selectedPoints;
        for (size_t i = 0; i < 6; i++) {
            for (size_t j = i + 1; j < 6; j++) {
                const T d = (vertexData_[extremeValues_[i]] - vertexData_[extremeValues_[j]]).squaredNorm();
                if (d > maxD) {
                    maxD = d;
                    selectedPoints = {extremeValues_[i], extremeValues_[j]};
                }
            }
        }
        if (maxD == epsilonSquared_) {
            // A degenerate case: the point cloud seems to consists of a single point
            return MeshBuilder<T>(0, std::min((size_t) 1, vertexCount - 1), std::min((size_t) 2, vertexCount - 1),
                                  std::min((size_t) 3, vertexCount - 1));
        }
        assert(selectedPoints.first != selectedPoints.second);

        // Find the most distant point to the line between the two chosen extreme points.
        const Ray<T> r(vertexData_[selectedPoints.first],
                       (vertexData_[selectedPoints.second] - vertexData_[selectedPoints.first]));
        maxD = epsilonSquared_;
        size_t maxI = std::numeric_limits<size_t>::max();
        const size_t vCount = vertexData_.size();
        for (size_t i = 0; i < vCount; i++) {
            const T distToRay = mathutils::getSquaredDistanceBetweenPointAndRay(vertexData_[i], r);
            if (distToRay > maxD) {
                maxD = distToRay;
                maxI = i;
            }
        }
        if (maxD == epsilonSquared_) {
            // It appears that the point cloud belongs to a 1 dimensional subspace of R^3: convex hull has no volume => return a thin triangle
            // Pick any point other than selectedPoints.first and selectedPoints.second as the third point of the triangle
            auto it = std::find_if(vertexData_.begin(), vertexData_.end(), [&](const Eigen::Matrix<T, 3, 1> &ve) {
                return ve != vertexData_[selectedPoints.first] && ve != vertexData_[selectedPoints.second];
            });
            const IndexType thirdPoint = (it == vertexData_.end()) ? selectedPoints.first : std::distance(
                    vertexData_.begin(), it);
            it = std::find_if(vertexData_.begin(), vertexData_.end(), [&](const Eigen::Matrix<T, 3, 1> &ve) {
                return ve != vertexData_[selectedPoints.first] && ve != vertexData_[selectedPoints.second] &&
                       ve != vertexData_[thirdPoint];
            });
            const IndexType fourthPoint = (it == vertexData_.end()) ? selectedPoints.first : std::distance(
                    vertexData_.begin(), it);
            return MeshBuilder<T>(selectedPoints.first, selectedPoints.second, thirdPoint, fourthPoint);
        }

        // These three points form the base triangle for our tetrahedron.
        assert(selectedPoints.first != maxI && selectedPoints.second != maxI);
        std::array<size_t, 3> baseTriangle{selectedPoints.first, selectedPoints.second, maxI};
        const Eigen::Matrix<T, 3, 1> baseTriangleVertices[] = {vertexData_[baseTriangle[0]],
                                                               vertexData_[baseTriangle[1]],
                                                               vertexData_[baseTriangle[2]]};

        // Next step is to find the 4th vertex of the tetrahedron. We naturally choose the point farthest away from the triangle plane.
        maxD = epsilon_;
        maxI = 0;
        const Eigen::Matrix<T, 3, 1> N = mathutils::getTriangleNormal(baseTriangleVertices[0], baseTriangleVertices[1],
                                                                      baseTriangleVertices[2]);
        Plane<T> trianglePlane(N, baseTriangleVertices[0]);
        for (size_t i = 0; i < vCount; i++) {
            const T d = std::abs(mathutils::getSignedDistanceToPlane(vertexData_[i], trianglePlane));
            if (d > maxD) {
                maxD = d;
                maxI = i;
            }
        }
        if (maxD == epsilon_) {
            // All the points seem to lie on a 2D subspace of R^3. How to handle this? Well, let's add one extra point to the point cloud so that the convex hull will have volume.
            planar_ = true;
            const Eigen::Matrix<T, 3, 1> N = mathutils::getTriangleNormal(baseTriangleVertices[1],
                                                                          baseTriangleVertices[2],
                                                                          baseTriangleVertices[0]);
            planarPointCloudTemp_.clear();
            planarPointCloudTemp_.insert(planarPointCloudTemp_.begin(), vertexData_.begin(), vertexData_.end());
            const Eigen::Matrix<T, 3, 1> extraPoint = N + vertexData_[0];
            planarPointCloudTemp_.push_back(extraPoint);
            maxI = planarPointCloudTemp_.size() - 1;
            vertexData_ = VertexDataSource<T>(planarPointCloudTemp_);
        }

        // Enforce CCW orientation (if user prefers clockwise orientation, swap two vertices in each triangle when final mesh is created)
        const Plane<T> triPlane(N, baseTriangleVertices[0]);
        if (triPlane.isPointOnPositiveSide(vertexData_[maxI])) {
            std::swap(baseTriangle[0], baseTriangle[1]);
        }

        // Create a tetrahedron half edge mesh and compute planes defined by each triangle
        MeshBuilder<T> mesh(baseTriangle[0], baseTriangle[1], baseTriangle[2], maxI);
        for (auto &f : mesh.faces_) {
            auto v = mesh.getVertexIndicesOfFace(f);
            const Eigen::Matrix<T, 3, 1> &va = vertexData_[v[0]];
            const Eigen::Matrix<T, 3, 1> &vb = vertexData_[v[1]];
            const Eigen::Matrix<T, 3, 1> &vc = vertexData_[v[2]];
            const Eigen::Matrix<T, 3, 1> N = mathutils::getTriangleNormal(va, vb, vc);
            const Plane<T> trianglePlane(N, va);
            f.P_ = trianglePlane;
        }

        // Finally we assign a face for each vertex outside the tetrahedron (vertices inside the tetrahedron have no role anymore)
        for (size_t i = 0; i < vCount; i++) {
            for (auto &face : mesh.faces_) {
                if (addPointToFace(face, i)) {
                    break;
                }
            }
        }
        return mesh;
    }

    /*
     * Explicit template specifications for float and double
     */

    template
    class QuickHull<float>;

    template
    class QuickHull<double>;
}

