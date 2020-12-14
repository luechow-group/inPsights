// SPDX-License-Identifier: CC0-1.0
// from https://github.com/akuukka/quickhull

#ifndef VERTEXDATASOURCE_H
#define VERTEXDATASOURCE_H

#include <Eigen/Core>

namespace quickhull {
	
	template<typename T>
	class VertexDataSource {
		const Eigen::Matrix<T,3,1>* ptr_;
		size_t count_;
	
	public:
		VertexDataSource(const Eigen::Matrix<T,3,1>* ptr, size_t count)
		: ptr_(ptr), count_(count) {
			
		}
		
		VertexDataSource(const std::vector<Eigen::Matrix<T,3,1>>& vec)
		: ptr_(&vec[0]), count_(vec.size()) {
			
		}
		
		VertexDataSource() : ptr_(nullptr), count_(0) {
			
		}
		
		VertexDataSource& operator=(const VertexDataSource& other) = default;
		
		size_t size() const {
			return count_;
		}
		
		const Eigen::Matrix<T,3,1>& operator[](size_t index) const {
			return ptr_[index];
		}
		
		const Eigen::Matrix<T,3,1>* begin() const {
			return ptr_;
		}
		
		const Eigen::Matrix<T,3,1>* end() const {
			return ptr_ + count_;
		}
	};
	
}

#endif //VERTEXDATASOURCE_H
