#ifndef Pool_h
#define Pool_h

#include <vector>
#include <memory>

namespace quickhull {
	
	template<typename T>
	class Pool {
		std::vector<std::unique_ptr<T>> data_;
	public:
		void clear() {
			data_.clear();
		}
		
		void reclaim(std::unique_ptr<T>& ptr) {
			data_.push_back(std::move(ptr));
		}
		
		std::unique_ptr<T> get() {
			if (data_.size()==0) {
				return std::unique_ptr<T>(new T());
			}
			auto it = data_.end()-1;
			std::unique_ptr<T> r = std::move(*it);
			data_.erase(it);
			return r;
		}
		
	};
}

#endif /* Pool_h */
