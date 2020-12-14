// SPDX-License-Identifier: CC0-1.0
// from https://github.com/akuukka/quickhull

#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include <Eigen/Core>
#include "Mesh.h"
#include "VertexDataSource.h"
#include <vector>
#include <unordered_map>
#include <fstream>
#include <memory>
#include <Vertex.h>
#include <Triangle.h>
#include <Enumerate.h>

namespace quickhull {

	template<typename T>
	class ConvexHull {
		std::unique_ptr<std::vector<Eigen::Matrix<T,3,1>>> optimizedVertexBuffer_;
		VertexDataSource<T> vertices_;
		std::vector<size_t> indices_;
	public:
		ConvexHull() = default;
		
		// Copy constructor
		ConvexHull(const ConvexHull& o) {
			indices_ = o.indices_;
			if (o.optimizedVertexBuffer_) {
				optimizedVertexBuffer_.reset(new std::vector<Eigen::Matrix<T,3,1>>(*o.optimizedVertexBuffer_));
				vertices_ = VertexDataSource<T>(*optimizedVertexBuffer_);
			}
			else {
				vertices_ = o.vertices_;
			}
		}
		
		ConvexHull& operator=(const ConvexHull& o) {
			if (&o == this) {
				return *this;
			}
			indices_ = o.indices_;
			if (o.optimizedVertexBuffer_) {
				optimizedVertexBuffer_.reset(new std::vector<Eigen::Matrix<T,3,1>>(*o.optimizedVertexBuffer_));
				vertices_ = VertexDataSource<T>(*optimizedVertexBuffer_);
			}
			else {
				vertices_ = o.vertices_;
			}
			return *this;
		}
		
		ConvexHull(ConvexHull&& o) {
			indices_ = std::move(o.indices_);
			if (o.optimizedVertexBuffer_) {
				optimizedVertexBuffer_ = std::move(o.optimizedVertexBuffer_);
				o.vertices_ = VertexDataSource<T>();
				vertices_ = VertexDataSource<T>(*optimizedVertexBuffer_);
			}
			else {
				vertices_ = o.vertices_;
			}
		}
		
		ConvexHull& operator=(ConvexHull&& o) {
			if (&o == this) {
				return *this;
			}
			indices_ = std::move(o.indices_);
			if (o.optimizedVertexBuffer_) {
				optimizedVertexBuffer_ = std::move(o.optimizedVertexBuffer_);
				o.vertices_ = VertexDataSource<T>();
				vertices_ = VertexDataSource<T>(*optimizedVertexBuffer_);
			}
			else {
				vertices_ = o.vertices_;
			}
			return *this;
		}
		
		// Construct vertex and index buffers from half edge mesh and pointcloud
		ConvexHull(const MeshBuilder<T>& mesh, const VertexDataSource<T>& pointCloud, bool CCW, bool useOriginalIndices) {
			if (!useOriginalIndices) {
				optimizedVertexBuffer_.reset(new std::vector<Eigen::Matrix<T,3,1>>());
			}
			
			std::vector<bool> faceProcessed(mesh.faces_.size(),false);
			std::vector<size_t> faceStack;
			std::unordered_map<size_t,size_t> vertexIndexMapping; // Map vertex indices from original point cloud to the new mesh vertex indices
			for (size_t i = 0;i<mesh.faces_.size();i++) {
				if (!mesh.faces_[i].isDisabled()) {
					faceStack.push_back(i);
					break;
				}
			}
			if (faceStack.empty()) {
				return;
			}
			
			const size_t finalMeshFaceCount = mesh.faces_.size() - mesh.disabledFaces_.size();
			indices_.reserve(finalMeshFaceCount*3);

			while (!faceStack.empty()) {
				auto it = faceStack.end()-1;
				size_t top = *it;
				assert(!mesh.faces_[top].isDisabled());
				faceStack.erase(it);
				if (faceProcessed[top]) {
					continue;
				}
				else {
					faceProcessed[top]=true;
					auto halfEdges = mesh.getHalfEdgeIndicesOfFace(mesh.faces_[top]);
					size_t adjacent[] = {mesh.halfEdges_[mesh.halfEdges_[halfEdges[0]].opp_].face_,mesh.halfEdges_[mesh.halfEdges_[halfEdges[1]].opp_].face_,mesh.halfEdges_[mesh.halfEdges_[halfEdges[2]].opp_].face_};
					for (auto a : adjacent) {
						if (!faceProcessed[a] && !mesh.faces_[a].isDisabled()) {
							faceStack.push_back(a);
						}
					}
					auto vertices = mesh.getVertexIndicesOfFace(mesh.faces_[top]);
					if (!useOriginalIndices) {
						for (auto& v : vertices) {
							auto foundIt = vertexIndexMapping.find(v);
							if (foundIt == vertexIndexMapping.end()) {
								optimizedVertexBuffer_->push_back(pointCloud[v]);
								vertexIndexMapping[v] = optimizedVertexBuffer_->size()-1;
								v = optimizedVertexBuffer_->size()-1;
							}
							else {
								v = foundIt->second;
							}
						}
					}
					indices_.push_back(vertices[0]);
					if (CCW) {
						indices_.push_back(vertices[2]);
						indices_.push_back(vertices[1]);
					}
					else {
						indices_.push_back(vertices[1]);
						indices_.push_back(vertices[2]);
					}
				}
			}
			
			if (!useOriginalIndices) {
				vertices_ = VertexDataSource<T>(*optimizedVertexBuffer_);
			}
			else {
				vertices_ = pointCloud;
			}
		}

		std::vector<size_t>& getIndexBuffer() {
            return indices_;
        }

        const std::vector<size_t>& getIndexBuffer() const {
            return indices_;
        }

		VertexDataSource<T>& getVertexBuffer() {
            return vertices_;
        }

        const VertexDataSource<T>& getVertexBuffer() const {
            return vertices_;
        }
		
		// Export the mesh to a Waveform OBJ file
		void writeWaveformOBJ(const std::string& filename, const std::string& objectName = "quickhull")
		{
			std::ofstream objFile;
			objFile.open (filename);
			objFile << "o " << objectName << "\n";
			for (const auto& v : getVertexBuffer()) {
				objFile << "v " << v.x() << " " << v.y() << " " << v.z() << "\n";
			}
			const auto& indBuf = getIndexBuffer();
			size_t triangleCount = indBuf.size()/3;
			for (size_t i=0;i<triangleCount;i++) {
				objFile << "f " << indBuf[i*3]+1 << " " << indBuf[i*3+1]+1 << " " << indBuf[i*3+2]+1 << "\n";
			}
			objFile.close();
		}

		std::vector<Vertex> getVertices() const {
		    std::vector<Vertex> vertices(getVertexBuffer().size());

            for (const auto& [i,v] : enumerate(getVertexBuffer()))
                vertices[i] = Vertex(v);

            return vertices;
		}

        std::vector<Triangle> getTriangles() const {
            std::vector<Triangle> triangles(getIndexBuffer().size()/3);

            for (size_t i = 0; i < getIndexBuffer().size()/3; ++i)
                triangles[i] = Triangle({
                    static_cast<int>(getIndexBuffer()[i*3+0]),
                    static_cast<int>(getIndexBuffer()[i*3+1]),
                    static_cast<int>(getIndexBuffer()[i*3+2])});

            return triangles;
        }
	};

}

#endif //CONVEXHULL_H
