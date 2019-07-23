#ifndef QHPLANE_HPP_
#define QHPLANE_HPP_

#include <Eigen/Core>

namespace quickhull {

	template<typename T>
	class Plane {
	public:
		Eigen::Matrix<T,3,1> m_N;
		
		// Signed distance (if normal is of length 1) to the plane from origin
		T m_D;
		
		// Normal length squared
		T m_sqrNLength;

		bool isPointOnPositiveSide(const Eigen::Matrix<T,3,1>& Q) const {
			T d = m_N.dot(Q)+m_D;
            return d >= 0;
        }

		Plane() = default;

		// Construct a plane using normal N and any point P on the plane
		Plane(const Eigen::Matrix<T,3,1>& N, const Eigen::Matrix<T,3,1>& P)
		: m_N(N), m_D(-N.dot(P)), m_sqrNLength(m_N.squaredNorm()) {
			
		}
	};

}


#endif /* PLANE_HPP_ */
