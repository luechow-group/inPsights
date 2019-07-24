#ifndef MESH_HPP_
#define MESH_HPP_

#include <vector>
#include <Eigen/Core>
#include "Plane.h"
#include "Pool.h"
#include "Types.h"
#include <string>
#include <array>
#include <cassert>
#include <limits>
#include <memory>
#include "VertexDataSource.h"
#include <unordered_map>

namespace quickhull {

	template <typename T>
	class MeshBuilder {
	public:
		struct HalfEdge {
			IndexType endVertex_;
			IndexType opp_;
			IndexType face_;
			IndexType next_;
			
			void disable() {
				endVertex_ = std::numeric_limits<IndexType>::max();
			}
			
			bool isDisabled() const {
				return endVertex_ == std::numeric_limits<IndexType>::max();
			}
		};

		struct Face {
			IndexType he_;
			Plane<T> P_;
			T mostDistantPointDist_;
			IndexType mostDistantPoint_;
			size_t visibilityCheckedOnIteration_;
			std::uint8_t isVisibleFaceOnCurrentIteration_ : 1;
			std::uint8_t inFaceStack_ : 1;
			std::uint8_t horizonEdgesOnCurrentIteration_ : 3; // Bit for each half edge assigned to this face, each being 0 or 1 depending on whether the edge belongs to horizon edge
			std::unique_ptr<std::vector<IndexType>> pointsOnPositiveSide_;

			Face() : he_(std::numeric_limits<IndexType>::max()),
					 mostDistantPointDist_(0),
					 mostDistantPoint_(0),
					 visibilityCheckedOnIteration_(0),
					 isVisibleFaceOnCurrentIteration_(0),
					 inFaceStack_(0),
					 horizonEdgesOnCurrentIteration_(0)
			{

			}
			
			void disable() {
				he_ = std::numeric_limits<IndexType>::max();
			}

			bool isDisabled() const {
				return he_ == std::numeric_limits<IndexType>::max();
			}
		};

		// Mesh data
		std::vector<Face> faces_;
		std::vector<HalfEdge> halfEdges_;
		
		// When the mesh is modified and faces and half edges are removed from it, we do not actually remove them from the container vectors.
		// Insted, they are marked as disabled which means that the indices can be reused when we need to add new faces and half edges to the mesh.
		// We store the free indices in the following vectors.
		std::vector<IndexType> disabledFaces_,disabledHalfEdges_;
		
		IndexType addFace() {
			if (disabledFaces_.size()) {
				IndexType index = disabledFaces_.back();
				auto& f = faces_[index];
				assert(f.isDisabled());
				assert(!f.pointsOnPositiveSide_);
				f.mostDistantPointDist_ = 0;
				disabledFaces_.pop_back();
				return index;
			}
			faces_.emplace_back();
			return faces_.size()-1;
		}

		IndexType addHalfEdge()	{
			if (disabledHalfEdges_.size()) {
				const IndexType index = disabledHalfEdges_.back();
				disabledHalfEdges_.pop_back();
				return index;
			}
			halfEdges_.emplace_back();
			return halfEdges_.size()-1;
		}

		// Mark a face as disabled and return a pointer to the points that were on the positive of it.
		std::unique_ptr<std::vector<IndexType>> disableFace(IndexType faceIndex) {
			auto& f = faces_[faceIndex];
			f.disable();
			disabledFaces_.push_back(faceIndex);
			return std::move(f.pointsOnPositiveSide_);
		}

		void disableHalfEdge(IndexType heIndex) {
			auto& he = halfEdges_[heIndex];
			he.disable();
			disabledHalfEdges_.push_back(heIndex);
		}

		MeshBuilder() = default;
		
		// Create a mesh with initial tetrahedron ABCD. Dot product of AB with the normal of triangle ABC should be negative.
		MeshBuilder(IndexType a, IndexType b, IndexType c, IndexType d) {
			// Create halfedges
			HalfEdge AB;
			AB.endVertex_ = b;
			AB.opp_ = 6;
			AB.face_ = 0;
			AB.next_ = 1;
			halfEdges_.push_back(AB);

			HalfEdge BC;
			BC.endVertex_ = c;
			BC.opp_ = 9;
			BC.face_ = 0;
			BC.next_ = 2;
			halfEdges_.push_back(BC);

			HalfEdge CA;
			CA.endVertex_ = a;
			CA.opp_ = 3;
			CA.face_ = 0;
			CA.next_ = 0;
			halfEdges_.push_back(CA);

			HalfEdge AC;
			AC.endVertex_ = c;
			AC.opp_ = 2;
			AC.face_ = 1;
			AC.next_ = 4;
			halfEdges_.push_back(AC);

			HalfEdge CD;
			CD.endVertex_ = d;
			CD.opp_ = 11;
			CD.face_ = 1;
			CD.next_ = 5;
			halfEdges_.push_back(CD);

			HalfEdge DA;
			DA.endVertex_ = a;
			DA.opp_ = 7;
			DA.face_ = 1;
			DA.next_ = 3;
			halfEdges_.push_back(DA);

			HalfEdge BA;
			BA.endVertex_ = a;
			BA.opp_ = 0;
			BA.face_ = 2;
			BA.next_ = 7;
			halfEdges_.push_back(BA);

			HalfEdge AD;
			AD.endVertex_ = d;
			AD.opp_ = 5;
			AD.face_ = 2;
			AD.next_ = 8;
			halfEdges_.push_back(AD);

			HalfEdge DB;
			DB.endVertex_ = b;
			DB.opp_ = 10;
			DB.face_ = 2;
			DB.next_ = 6;
			halfEdges_.push_back(DB);

			HalfEdge CB;
			CB.endVertex_ = b;
			CB.opp_ = 1;
			CB.face_ = 3;
			CB.next_ = 10;
			halfEdges_.push_back(CB);

			HalfEdge BD;
			BD.endVertex_ = d;
			BD.opp_ = 8;
			BD.face_ = 3;
			BD.next_ = 11;
			halfEdges_.push_back(BD);

			HalfEdge DC;
			DC.endVertex_ = c;
			DC.opp_ = 4;
			DC.face_ = 3;
			DC.next_ = 9;
			halfEdges_.push_back(DC);

			// Create faces
			Face ABC;
			ABC.he_ = 0;
			faces_.push_back(std::move(ABC));

			Face ACD;
			ACD.he_ = 3;
			faces_.push_back(std::move(ACD));

			Face BAD;
			BAD.he_ = 6;
			faces_.push_back(std::move(BAD));

			Face CBD;
			CBD.he_ = 9;
			faces_.push_back(std::move(CBD));
		}

		std::array<IndexType,3> getVertexIndicesOfFace(const Face& f) const {
			std::array<IndexType,3> v;
			const HalfEdge* he = &halfEdges_[f.he_];
			v[0] = he->endVertex_;
			he = &halfEdges_[he->next_];
			v[1] = he->endVertex_;
			he = &halfEdges_[he->next_];
			v[2] = he->endVertex_;
			return v;
		}

		std::array<IndexType,2> getVertexIndicesOfHalfEdge(const HalfEdge& he) const {
			return {halfEdges_[he.opp_].endVertex_,he.endVertex_};
		}

		std::array<IndexType,3> getHalfEdgeIndicesOfFace(const Face& f) const {
			return {f.he_,halfEdges_[f.he_].next_,halfEdges_[halfEdges_[f.he_].next_].next_};
		}
	};
}

#endif 
