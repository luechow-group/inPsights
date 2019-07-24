
#ifndef QuickHull_MathUtils_hpp
#define QuickHull_MathUtils_hpp

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Ray.h"

namespace quickhull {
	
	namespace mathutils {
		
		template <typename T>
		inline T getSquaredDistanceBetweenPointAndRay(const Eigen::Matrix<T,3,1>& p, const Ray<T>& r) {
			const Eigen::Matrix<T,Eigen::Dynamic,1> s = p-r.S_;
			T t = s.dot(r.V_);
			return s.squaredNorm() - t*t*r.VInvLengthSquared_;
		}
		
		// Note that the unit of distance returned is relative to plane's normal's length (divide by N.getNormalized() if needed to get the "real" distance).
		template <typename T>
		inline T getSignedDistanceToPlane(const Eigen::Matrix<T,3,1>& v, const Plane<T>& p) {
			return p.N_.dot(v) + p.D_;
		}

		template <typename T>
		inline Eigen::Matrix<T,3,1> getTriangleNormal(
		        const Eigen::Matrix<T,3,1>& a,
		        const Eigen::Matrix<T,3,1>& b,
		        const Eigen::Matrix<T,3,1>& c) {
		    return (a-c).cross(b-c);
		}
	}
	
}


#endif
