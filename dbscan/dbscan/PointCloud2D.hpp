#pragma once
#include <vector>
#include <stdexcept>
#include <set>

namespace dm {
	namespace dbscan {
		template <typename T>
		struct PointCloud2D{
			static constexpr int Noise = 0;
			static constexpr int Undefined = -1;
			struct Point{
				T  x, y;
				int label = Undefined;// cluster id
			};

			using ClusterSet = std::set<int>;

			std::vector<Point>  pts;
			int nClusters = 0;

			// Must return the number of data points
			inline size_t kdtree_get_point_count() const { return pts.size(); }
			inline size_t size()const { return pts.size(); }

			// Returns the dim'th component of the idx'th point in the class:
			// Since this is inlined and the "dim" argument is typically an immediate value, the
			//  "if/else's" are actually solved at compile time.
			inline T kdtree_get_pt(const size_t idx, const size_t dim) const
			{
				if (dim == 0) return pts[idx].x;
				else if (dim == 1) return pts[idx].y;
				else throw std::logic_error("Invalid dimension for 2d point cloud");
			}

			// Optional bounding-box computation: return false to default to a standard bbox computation loop.
			//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
			//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
			template <class BBOX>
			bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

			int label(const size_t idx)const { return pts[idx].label; }
			void setLabel(const size_t idx, int label) { pts[idx].label = label; }
			const Point& getPt( const size_t idx )const { return pts[idx]; }
			void setClusterNum(int n) { nClusters = n; }
			int getClusterNum()const { return nClusters; }

			ClusterSet getClusterIds() {
				ClusterSet clusters;
				for (size_t idx = 0; idx < size(); idx++) {
					if (pts[idx].label != Undefined)
						clusters.insert(pts[idx].label);
				}
				return clusters;
			}

			std::vector<Point> getPointsByClusterId(int cid) {
				std::vector<Point> points;
				for (size_t idx = 0; idx < pts.size(); idx++) {
					if (pts[idx].label == cid)
						points.push_back(pts[idx]);
				}
				return points;
			}
		};
	}// namespace dbscan
}// namespace dm
