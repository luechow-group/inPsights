#ifndef QuickHull_Ray_hpp
#define QuickHull_Ray_hpp

#include <Eigen/Core>

namespace quickhull {

	template <typename T>
	struct Ray {
		const Eigen::Matrix<T,3,1> m_S;
		const Eigen::Matrix<T,3,1> m_V;
		const T m_VInvLengthSquared;
		
		Ray(const Eigen::Matrix<T,3,1>& S,const Eigen::Matrix<T,3,1>& V)
		: m_S(S), m_V(V), m_VInvLengthSquared(1/m_V.squaredNorm()) {
		}
	};
	
}


#endif
