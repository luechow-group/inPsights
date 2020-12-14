// SPDX-License-Identifier: CC0-1.0
// from https://github.com/akuukka/quickhull

#ifndef RAY_H
#define RAY_H

#include <Eigen/Core>

namespace quickhull {

	template <typename T>
	struct Ray {
		const Eigen::Matrix<T,3,1> S_;
		const Eigen::Matrix<T,3,1> V_;
		const T VInvLengthSquared_;
		
		Ray(const Eigen::Matrix<T,3,1>& S,const Eigen::Matrix<T,3,1>& V)
		: S_(S), V_(V), VInvLengthSquared_(1/V_.squaredNorm()) {
		}
	};
	
}

#endif //RAY_H
