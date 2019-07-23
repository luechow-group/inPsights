#ifndef VertexDataSource_h
#define VertexDataSource_h

#include <Eigen/Core>

namespace quickhull {
	
	template<typename T>
	class VertexDataSource {
		const Eigen::Matrix<T,3,1>* m_ptr;
		size_t m_count;
	
	public:
		VertexDataSource(const Eigen::Matrix<T,3,1>* ptr, size_t count)
		: m_ptr(ptr), m_count(count) {
			
		}
		
		VertexDataSource(const std::vector<Eigen::Matrix<T,3,1>>& vec)
		: m_ptr(&vec[0]), m_count(vec.size()) {
			
		}
		
		VertexDataSource() : m_ptr(nullptr), m_count(0) {
			
		}
		
		VertexDataSource& operator=(const VertexDataSource& other) = default;
		
		size_t size() const {
			return m_count;
		}
		
		const Eigen::Matrix<T,3,1>& operator[](size_t index) const {
			return m_ptr[index];
		}
		
		const Eigen::Matrix<T,3,1>* begin() const {
			return m_ptr;
		}
		
		const Eigen::Matrix<T,3,1>* end() const {
			return m_ptr + m_count;
		}
	};
	
}


#endif /* VertexDataSource_h */
