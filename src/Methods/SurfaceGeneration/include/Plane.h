#ifndef QHPLANE_HPP_
#define QHPLANE_HPP_

#include <Eigen/Core>

namespace quickhull {

	template<typename T>
	class Plane {
	public:
		Eigen::Matrix<T,3,1> N_;
		
		// Signed distance (if normal is of length 1) to the plane from origin
		T D_;
		
		// Normal length squared
		T sqrNLength_;

		bool isPointOnPositiveSide(const Eigen::Matrix<T,3,1>& Q) const {
			T d = N_.dot(Q)+D_;
            return d >= 0;
        }

		Plane() = default;

		// Construct a plane using normal N and any point P on the plane
		Plane(const Eigen::Matrix<T,3,1>& N, const Eigen::Matrix<T,3,1>& P)
		: N_(N), D_(-N.dot(P)), sqrNLength_(N_.squaredNorm()) {
			
		}
	};

}


#endif /* PLANE_HPP_ */
