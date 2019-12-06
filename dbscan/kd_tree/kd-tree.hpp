/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
 * Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
 * Copyright 2011-2016  Jose Luis Blanco (joseluisblancoc@gmail.com).
 *   All rights reserved.
 *
 * THE BSD LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

 /** \mainpage nanoflann C++ API documentation
  *  nanoflann is a C++ header-only library for building KD-Trees, mostly
  *  optimized for 2D or 3D point clouds.
  *
  *  nanoflann does not require compiling or installing, just an
  *  #include <nanoflann.hpp> in your code.
  *
  *  See:
  *   - <a href="modules.html" >C++ API organized by modules</a>
  *   - <a href="https://github.com/jlblancoc/nanoflann" >Online README</a>
  *   - <a href="http://jlblancoc.github.io/nanoflann/" >Doxygen
  * documentation</a>
  */
#pragma once
#include "PoolAllocator.hpp"

#include <array>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <cassert>

namespace dm {
	namespace KdTree {
		/** kd-tree base-class
		 *
		 * Contains the member functions common to the classes KDTreeSingleIndexAdaptor
		 * and KDTreeSingleIndexDynamicAdaptor_.
		 *
		 * \tparam Derived The name of the class which inherits this class.
		 * \tparam DatasetAdaptor The user-provided adaptor (see comments above).
		 * \tparam Distance The distance metric to use, these are all classes derived
		 * from nanoflann::Metric \tparam DIM Dimensionality of data points (e.g. 3 for
		 * 3D points) \tparam IndexType Will be typically size_t or int
		*/
		template<
			typename Derived, 
			typename Distance, 
			typename DatasetAdaptor, 
			int DIM, 
			typename IndexType = size_t
		>
		class KDTreeBaseClass {
		public: // == TYPES ==
			typedef typename Distance::ElementType ElementType;
			typedef typename Distance::DistanceType DistanceType;
			
			struct Node {
				/** Union used because a node can be either a LEAF node or a non-leaf node,
				  * so both data fields are never used simultaneously 
				  */
				union {
					struct leaf{
						IndexType left, right;//!< Indices of points in leaf node
					} lr;
					struct nonleaf {
						int divfeat;					//!< Dimension used for subdivision.
						DistanceType divlow, divhigh;	//!< The values used for subdivision.
					} sub;
				} node_type;
				Node* child1;
				Node* child2;
			};

			typedef Node* NodePtr;
			struct Interval {
				ElementType low, high;
			};

			template <int DIM, typename T> struct array_or_vector_selector {
				typedef std::array<T, DIM> container_t;
			};
			/** Dynamic size version */
			template <typename T> struct array_or_vector_selector<-1, T> {
				typedef std::vector<T> container_t;
			};

			/** Define "BoundingBox" as a fixed-size or variable-size container depending
			 *  on "DIM" 
			 */
			typedef typename array_or_vector_selector<DIM, Interval>::container_t BoundingBox;
			/** Define "distance_vector_t" as a fixed-size or variable-size container
			  * depending on "DIM" 
			  */
			typedef typename array_or_vector_selector<DIM, DistanceType>::container_t distance_vector_t;
		public: // == CONSTANTS ==
			static constexpr DistanceType Eps = 0.00001;
		public: // == PUBLIC MEMBERS ==
			std::vector<IndexType>	vind;
			NodePtr					root_node;
			size_t					m_leaf_max_size;
			size_t					m_size;				//!< Number of current points in the dataset
			size_t					m_size_at_index_build;	//!< Number of points in the dataset when the
														//!< index was built
			int						dim;				//!< Dimensionality of each data point
			BoundingBox				root_bbox;			//!< The KD-tree used to find neighbours
			/**
			* Pooled memory allocator.
			*
			* Using a pooled memory allocator is more efficient
			* than allocating memory directly when there is a large
			* number small of memory allocations.
			*/
			PoolAllocator			pool;
		public:	// == FUNCTIONS ==
			 /** Returns number of points in dataset  */
			size_t size(const Derived &obj) const { return obj.m_size; }

			/** Returns the length of each point in the dataset */
			size_t veclen(const Derived &obj) {
				return static_cast<size_t>(DIM > 0 ? DIM : obj.dim);
			}

			/** Frees the previously-built index. Automatically called within
			  * buildIndex(). 
			  */
			void freeIndex(Derived &obj) {
				obj.pool.freeAll();
				obj.root_node = nullptr;
				obj.m_size_at_index_build = 0;
			}

			inline ElementType dataset_get(const Derived& obj, size_t idx, int dim)const {
				return obj.dataset.kdtree_get_pt(idx, dim);
			}

			void computeMinMax(const Derived& obj, IndexType* ind, IndexType count, int dim, ElementType& minElem, ElementType& maxElem) {
				minElem = dataset_get(obj, ind[0], dim);
				maxElem = dataset_get(obj, ind[0], dim);

				for (IndexType i = 1; i < count; i++) {
					ElementType val = dataset_get(obj, ind[i], dim);

					if (val < minElem) minElem = val;

					if (val > maxElem) maxElem = val;
				}
			}

			/**
			 *  Subdivide the list of points by a plane perpendicular on axe corresponding
			 *  to the 'cutfeat' dimension at 'cutval' position.
			 *
			 *  On return:
			 *  dataset[ind[0..lim1-1]][cutfeat]<cutval
			 *  dataset[ind[lim1..lim2-1]][cutfeat]==cutval
			 *  dataset[ind[lim2..count]][cutfeat]>cutval
			 */
			void planeSplit(Derived &obj, IndexType *ind, const IndexType count, int cutfeat, DistanceType &cutval, IndexType &lim1, IndexType &lim2) {
				IndexType left = 0;
				IndexType right = count - 1;
				
				for (;; ) {
					while (left <= right && dataset_get(obj, ind[left], cutfeat) < cutval) {
						++left;
					}
					while (right && left <= right && dataset_get(obj, ind[right], cutfeat) >= cutval) {
						--right;
					}

					if (left > right || !right) {
						break;
					}

					std::swap(ind[left], ind[right]);
					++left;
					--right;
				}
				/* If either list is empty, it means that all remaining features
				 * are identical. Split in the middle to maintain a balanced tree.
				 */
				lim1 = left;
				right = count - 1;
				for (;;) {
					while (left <= right && dataset_get(obj, ind[left], cutfeat) <= cutval)
						++left;
					while (right && left <= right &&
						dataset_get(obj, ind[right], cutfeat) > cutval)
						--right;
					if (left > right || !right)
						break; // "!right" was added to support unsigned Index types
					std::swap(ind[left], ind[right]);
					++left;
					--right;
				}
				lim2 = left;
			}

			/**
			  * Create a tree node that subdivides the list of vecs from vind[first]
			  * to vind[last].  The routine is called recursively on each sublist.
			  *
			  * @param left index of the first vector
			  * @param right index of the last vector
			  */
			void middleSplit(Derived &obj, IndexType *ind, IndexType count, IndexType &index, int &cutfeat, DistanceType &cutval, const BoundingBox &bbox) {
				ElementType maxspan = bbox[0].high - bbox[0].low;

				for (int idx = 1; idx < DIM; idx++) {
					ElementType span = bbox[idx].high - bbox[idx].low;

					if (span > maxspan) maxspan = span;
				}

				ElementType maxSpread = -1;
				cutfeat = 0;

				for (int idx = 0; idx < DIM; idx++) {
					ElementType span = bbox[idx].high - bbox[idx].low;

					if (span > (1 - Eps) * maxspan) {
						volatile bool stope = true;
						ElementType minE, maxE;
						computeMinMax(obj, ind, count, idx, minE, maxE);
						ElementType spread = maxE - minE;
						
						if (spread > maxSpread) {
							cutfeat = idx;
							maxSpread = spread;
						}
					}
				}

				// split in the middle
				DistanceType splitVal = (bbox[cutfeat].low + bbox[cutfeat].high) / 2;
				ElementType minElem, maxElem;
				computeMinMax(obj, ind, count, cutfeat, minElem, maxElem);

				if (splitVal < minElem)
					cutval = minElem;
				else if (splitVal > maxElem)
					cutval = maxElem;
				else
					cutval = splitVal;

				IndexType lim1, lim2;
				planeSplit(obj, ind, count, cutfeat, cutval, lim1, lim2);

				if (lim1 > count / 2)
					index = lim1;
				else if (lim2 < count / 2)
					index = lim2;
				else
					index = count / 2;
			}

			/**
			 * Create a tree node that subdivides the list of vecs from vind[first]
			 * to vind[last].  The routine is called recursively on each sublist.
			 *
			 * @param left index of the first vector
			 * @param right index of the last vector
			 */
			NodePtr divideTree(Derived &obj, const IndexType left, const IndexType right, BoundingBox &bbox) {
				NodePtr node = obj.pool.template allocate<Node>();

				if ((right - left) <= static_cast<IndexType>(obj.m_leaf_max_size)) {
					node->child1 = node->child2 = nullptr;// mark as leaf node
					node->node_type.lr.left = left;
					node->node_type.lr.right = right;

					// compute bounding box of leaf point
					for (int idx = 0; idx < DIM; idx++) {
						bbox[idx].low = dataset_get(obj, obj.vind[left], idx);
						bbox[idx].high = dataset_get(obj, obj.vind[left], idx);
					}

					for (IndexType k = left + 1; k < right; k++) {
						for (int d = 0; d < DIM; d++) {
							if (bbox[d].low > dataset_get(obj, obj.vind[k], d)) bbox[d].low = dataset_get(obj, obj.vind[k], d);
							if (bbox[d].high < dataset_get(obj, obj.vind[k], d)) bbox[d].high = dataset_get(obj, obj.vind[k], d);
						}
					}
				}else {
					IndexType idx;
					int cutfeat;
					DistanceType cutval;
					middleSplit(obj, &obj.vind[0] + left, right - left, idx, cutfeat, cutval, bbox);
					node->node_type.sub.divfeat = cutfeat;
					
					BoundingBox left_bbox(bbox);
					left_bbox[cutfeat].high = cutval;
					node->child1 = divideTree(obj, left, left + idx, left_bbox);

					BoundingBox right_box(bbox);
					right_box[cutfeat].low = cutval;
					node->child2 = divideTree(obj, left + idx, right, right_box);

					node->node_type.sub.divlow = left_bbox[cutfeat].high;// value, used for subdivision
					node->node_type.sub.divhigh = right_box[cutfeat].low;// value, used for subdivision
					
					for (int d = 0; d < DIM; d++) {
						bbox[d].low = std::min(left_bbox[d].low, right_box[d].low);
						bbox[d].high = std::max(left_bbox[d].high, right_box[d].high);
					}

					//volatile bool tmp = true;
				}

				return node;
			}
		
			/*
			 * count distances if point lays out of data area.
			*/
			DistanceType computeInitialDistances(const Derived& obj, const ElementType* vec, distance_vector_t& dists)const {
				assert(vec);
				DistanceType distsq = DistanceType();

				for (int d = 0; d < DIM; d++) {
					if (vec[d] < obj.root_bbox[d].low) {
						dists[d] = obj.distance.accum_dist(vec[d], obj.root_bbox[d].low, d);
						distsq += dists[d];
					}

					if (vec[d] > obj.root_bbox[d].high) {
						dists[d] = obj.distance.accum_dist(vec[d], obj.root_bbox[d].high, d);
						distsq += dists[d];
					}
				}

				return distsq;
			}
		};

		//-------------------------------------------------------------------------------------------
		//						KDTreeSingleIndexAdaptor
		//-------------------------------------------------------------------------------------------
		struct KDTreeSingleIndexAdaptorParams {
			KDTreeSingleIndexAdaptorParams(size_t _leaf_max_size = 10) : leaf_max_size(_leaf_max_size) {}
			size_t leaf_max_size;
		};

		/** Search options for KDTreeSingleIndexAdaptor::findNeighbors() */
		struct SearchParams {
			int checks = 0;		//!< Ignored parameter (Kept for compatibility with the FLANN
								//!< interface).
			float eps = 0.0f;	//!< search for eps-approximate neighbours (default: 0)
			bool sorted = true;	//!< only for radius search, require neighbours sorted by

								//!< distance (default: true)
			SearchParams(int checks_IGNORED_ = 32, float epsilon = 0, bool s = true)
				: checks( checks_IGNORED_ )
				, eps( epsilon )
				, sorted( s )
			{}
		};

		struct IndexDist_Sorter {
			/** PairType will be typically: std::pair<IndexType,DistanceType> */
			template <typename PairType>
			inline bool operator()(const PairType &p1, const PairType &p2) const {
				return p1.second < p2.second;
			}
		};

		/**
		 * A result-set class used when performing a radius based search.
		 */
		template<typename DT, typename IT = size_t>
		class RadiusResultSet {
		public: // == TYPES ==
			typedef DT DistanceType;
			typedef IT IndexType;
		public: // == MEMBERS ==
			const DistanceType										radius;
			std::vector< std::pair< IndexType, DistanceType > >&	m_indices_dists;
		public: // == CTOR ==
			RadiusResultSet(DistanceType r, std::vector< std::pair< IndexType, DistanceType > >& indices)
				: radius(r)
				, m_indices_dists(indices) {
				init();
			}

		public: // == METHODS ==
			void init() { clear(); }
			
			void clear() { m_indices_dists.clear(); }
			
			size_t size()const { return m_indices_dists.size(); }
			
			bool fool()const { return true; }
			
			bool addPoint(DistanceType dist, IndexType index) {
				if (dist < radius) {
					m_indices_dists.push_back({ index, dist });
				}

				return true;
			}
			
			DistanceType worstDist()const { return radius; }

			std::pair<IndexType, DistanceType> worst_item() const {
				if (m_indices_dists.empty())
					throw std::runtime_error("Can't invoke worst_item, indices empty");
				using PT = std::pair<IndexType, DistanceType>;
				auto it = std::max_element(m_indices_dists.begin(), m_indices_dists.end(), IndexDist_Sorter());
				return *it;
			}
		};

		template<typename T, typename = int> 
		struct has_resize : std::false_type{};

		template<typename T>
		struct has_resize<T, decltype( (void)std::declval<T>().resize(1), 0)>
			: std::true_type{};

		template<typename T, typename = int>
		struct has_assign : std::false_type {};

		template<typename T>
		struct has_assign<T, decltype((void)std::declval<T>().assign(1, 0), 0)> : std::true_type{};

		/**
		 * Free function to resize a resizable object
		 */
		template<typename Cnt>
		inline typename std::enable_if<has_resize<Cnt>::value, void>::type resize(Cnt& c, const size_t num) {
			c.resize(num);
		}

		/**
		 * Free function that has no effects on non resizable containers (e.g.
		 * std::array) It raises an exception if the expected size does not match
		 */
		template <typename Container>
		inline typename std::enable_if<!has_resize<Container>::value, void>::type resize(Container &c, const size_t nElements) {
			if (nElements != c.size())
				throw std::logic_error("Try to change the size of a std::array.");
		}

		/**
		 * Free function to assign to a container
		 */
		template<typename Container, typename T>
		inline typename std::enable_if<has_assign<Container>::value, void>::type assign(Container& c, const size_t nElements, const T& value) {
			c.assign(nElements, value);
		}

		/**
		 * Free function to assign to a std::array
		 */
		template <typename Container, typename T>
		inline typename std::enable_if<!has_assign<Container>::value, void>::type
		assign(Container &c, const size_t nElements, const T &value) {
			assert(c.size() >= nElements);
			for (size_t i = 0; i < nElements; i++)
				c[i] = value;
		}

		template<
			typename Distance,
			typename DatasetAdaptor,
			int DIM,
			typename IndexType = size_t
		>
		class KDTreeSingleIndexAdaptor : public KDTreeBaseClass
		<
			KDTreeSingleIndexAdaptor<Distance, DatasetAdaptor, DIM, IndexType>,
			Distance,
			DatasetAdaptor,
			DIM,
			IndexType
		>
		{
		public: // == TYPES ==
			typedef typename KDTreeBaseClass
			<
				KDTreeSingleIndexAdaptor
				<
					Distance, 
					DatasetAdaptor, 
					DIM,
					IndexType
				>,
				Distance, 
				DatasetAdaptor, 
				DIM, 
				IndexType
			> BaseClassRef;
			
			typedef typename BaseClassRef::ElementType ElementType;
			
			typedef typename BaseClassRef::DistanceType DistanceType;

			typedef typename BaseClassRef::Node Node;

			typedef Node *NodePtr;

			typedef typename BaseClassRef::Interval Interval;

			/** Define "BoundingBox" as a fixed-size or variable-size container depending
			 * on "DIM" 
			 */
			typedef typename BaseClassRef::BoundingBox BoundingBox;

			/** Define "distance_vector_t" as a fixed-size or variable-size container
			 * depending on "DIM" 
			 */
			typedef typename BaseClassRef::distance_vector_t distance_vector_t;

		public: // == MEMBERS ==
			const DatasetAdaptor&					dataset;
			const KDTreeSingleIndexAdaptorParams	index_params;
			Distance								distance;

		public: // == CTORs ==
			KDTreeSingleIndexAdaptor(const KDTreeSingleIndexAdaptor<Distance, DatasetAdaptor, DIM, IndexType>&) = delete;

		  /**
			* KDTree constructor
			*
			* Refer to docs in README.md or online in
			* https://github.com/jlblancoc/nanoflann
			*
			* The KD-Tree point dimension (the length of each point in the datase, e.g. 3
			* for 3D points) is determined by means of:
			*  - The \a DIM template parameter if >0 (highest priority)
			*  - Otherwise, the \a dimensionality parameter of this constructor.
			*
			* @param inputData Dataset with the input features
			* @param params Basically, the maximum leaf node size
			*/
			KDTreeSingleIndexAdaptor(const int dimensionality,
				const DatasetAdaptor &inputData,
				const KDTreeSingleIndexAdaptorParams &params =
				KDTreeSingleIndexAdaptorParams())
				: dataset(inputData)
				, index_params(params)
				, distance(inputData) {
				static_assert(DIM > 0, "Invalid dimension");
				BaseClassRef::root_node = nullptr;
				BaseClassRef::m_size = dataset.kdtree_get_point_count();
				BaseClassRef::m_size_at_index_build = BaseClassRef::m_size;
				BaseClassRef::dim = DIM;
				BaseClassRef::m_leaf_max_size = params.leaf_max_size;
				// Create a permutable array of indices to the input vectors.
				init_vind();
			}
		public: // == FUNCTIONS ==
			 /** Make sure the auxiliary list \a vind has the same size than the current
			  * dataset, and re-generate if size has changed. 
			  */
			void init_vind() {
				BaseClassRef::m_size = dataset.kdtree_get_point_count();

				if (BaseClassRef::vind.size() != BaseClassRef::m_size)
					BaseClassRef::vind.resize(BaseClassRef::m_size);

				for (size_t idx = 0; idx < BaseClassRef::m_size; idx++) {
					BaseClassRef::vind[idx] = idx;
				}
			}

			void computeBoundingBox(BoundingBox &bbox) {
				resize(bbox, DIM);

				if ( !dataset.kdtree_get_bbox(bbox)) {
					const size_t N = this->dataset.kdtree_get_point_count();

					if (!N)
						throw std::runtime_error("can't compute bounding box, no points specified");

					for (int d = 0; d < DIM; d++) {
						bbox[d].low = bbox[d].high = this->dataset_get(*this, 0, d);
					}

					for (size_t k = 1; k < N; k++) {
						for (int d = 0; d < DIM; d++) {
							if (this->dataset_get(*this, k, d) < bbox[d].low)
								bbox[d].low = this->dataset_get(*this, k, d);
							
							if (this->dataset_get(*this, k, d) > bbox[d].high)
								bbox[d].high = this->dataset_get(*this, k, d);
						}
					}
				}
			}

			/**
			 * Builds the index
			 */
			void buildIndex() {
				BaseClassRef::m_size = dataset.kdtree_get_point_count();
				BaseClassRef::m_size_at_index_build = BaseClassRef::m_size;
				init_vind();
				this->freeIndex( *this );
				BaseClassRef::m_size_at_index_build = BaseClassRef::m_size;
				
				if (BaseClassRef::m_size == 0)
					return;

				computeBoundingBox(BaseClassRef::root_bbox);
				BaseClassRef::root_node = this->divideTree(*this, 0, BaseClassRef::m_size, BaseClassRef::root_bbox);
			}
		
			/**
			* Performs an exact search in the tree starting from a node.
			* \tparam RESULTSET Should be any ResultSet<DistanceType>
			* \return true if the search should be continued, false if the results are
			* sufficient
			*/
			template<typename ResultSet>
			bool searchLevel(ResultSet& result_set, const ElementType* vec, const NodePtr node, DistanceType mindistq, distance_vector_t& dists, const float epsError)const {
				/* If this is a leaf node, then do check and return. */
				if (node->child1 == nullptr && node->child2 == nullptr) {
					DistanceType worst_dist = result_set.worstDist();

					for (IndexType i = node->node_type.lr.left; i < node->node_type.lr.right; ++i) {
						const IndexType index = BaseClassRef::vind[i];
						DistanceType dist = distance.evalMetric(vec, index, DIM);

						if (dist < worst_dist) {
							if (!result_set.addPoint(dist, BaseClassRef::vind[i])) {
								// the resultset doesn't want to receive any more points, we're done
								// searching!
								return false;
							}
						}
					}
					return true;
				}// isLeaf ENDed

				/* Which child branch should be taken first? */
				int idx = node->node_type.sub.divfeat;
				ElementType val = vec[idx];// value of query { 1, 1, 1 }, dimension
				DistanceType diff1 = val - node->node_type.sub.divlow;//val - lefbox[high]
				DistanceType diff2 = val - node->node_type.sub.divhigh;//val - right_box[low]
				NodePtr bestChild;
				NodePtr otherChild;
				DistanceType cut_dist;// distance to otherChild

				if (diff1 + diff2 < 0) {
					bestChild = node->child1;
					otherChild = node->child2;
					cut_dist = distance.accum_dist(val, node->node_type.sub.divhigh, idx);// right_box[low], dist to other child
				}else {
					bestChild = node->child2;
					otherChild = node->child1;
					cut_dist = distance.accum_dist(val, node->node_type.sub.divlow, idx);// lefbox[high]
				}
				
				/* Call recursively to search next level down. */
				if (!searchLevel(result_set, vec, bestChild, mindistq, dists, epsError) ) {
					return false;
				}

				DistanceType dst = dists[ idx ];
				// mindistq > 0 if query point lays out of area, sum of squared distances of all dimensions, euclidian dist
				// cut_dist - dist to otherChild
				// dists - distance squares from point to low, high, one per dimension. dist -> dim[d].low
				// euclidian_dist_to_other_child = euclidian dist + dist to other child  - distance_of_cur_dimension ^ 2
				mindistq = mindistq + cut_dist - dst;
				dists[idx] = cut_dist;

				if (mindistq * epsError <= result_set.worstDist()) {// should I check otherChild?
					
					if (!searchLevel(result_set, vec, otherChild, mindistq, dists, epsError)) {
						// the resultset doesn't want to receive any more points, we're done
						// searching!
						return false;
					}
				}
				dists[idx] = dst;
				return true;
			}// search level

			/** \name Query methods
			/**
			 * Find set of nearest neighbors to vec[0:dim-1]. Their indices are stored
			 * inside the result object.
			 *
			 * Params:
			 *     result = the result object in which the indices of the
			 * nearest-neighbors are stored vec = the vector for which to search the
			 * nearest neighbors
			 *
			 * \tparam RESULTSET Should be any ResultSet<DistanceType>
			 * \return  True if the requested neighbors could be found.
			 * \sa knnSearch, radiusSearch
			 */

			template<typename ResultSet>
			bool findNeighbors(ResultSet& result, const ElementType* vec, const SearchParams& searchParams)const {
				assert(vec);

				if (this->size(*this) == 0)
					return false;

				if (!BaseClassRef::root_node)
					throw std::runtime_error("findNeighbors called before building the index.");

				float epsError = 1 + searchParams.eps;
				distance_vector_t dists;
				auto zero = static_cast<decltype(result.worstDist())>(0);
				assign(dists, DIM, zero);
				DistanceType distsq = this->computeInitialDistances(*this, vec, dists);
				searchLevel(result, vec, BaseClassRef::root_node, distsq, dists, epsError);
				return result.fool();
			}

			/**
			 * Just like radiusSearch() but with a custom callback class for each point
			 * found in the radius of the query. See the source of RadiusResultSet<> as a
			 * start point for your own classes. \sa radiusSearch
			 */
			template<typename SearchCallback>
			size_t radiusSearchCustomCallback(const ElementType* query_point, SearchCallback & resultSet, const SearchParams& searchParams = SearchParams())const {
				this->findNeighbors(resultSet, query_point, searchParams);
				return resultSet.size();
			}

			/**
			 * Find all the neighbors to \a query_point[0:dim-1] within a maximum radius.
			 *  The output is given as a vector of pairs, of which the first element is a
			 *  point index and the second the corresponding distance. Previous contents of
			 * \a IndicesDists are cleared.
			 *
			 *  If searchParams.sorted==true, the output list is sorted by ascending
			 * distances.
			 *
			 *  For a better performance, it is advisable to do a .reserve() on the vector
			 * if you have any wild guess about the number of expected matches.
			 *
			 *  \sa knnSearch, findNeighbors, radiusSearchCustomCallback
			 * \return The number of points within the given radius (i.e. indices.size()
			 * or dists.size() )
			*/
			size_t radiusSearch(const ElementType* query_point, const DistanceType& radius, std::vector< std::pair< IndexType, DistanceType > >& IndicesDists, const SearchParams& searchParams)const {
				RadiusResultSet<DistanceType, IndexType> resultSet(radius, IndicesDists);
				const size_t nFound = radiusSearchCustomCallback(query_point, resultSet, searchParams);

				if (searchParams.sorted)
					std::sort(IndicesDists.begin(), IndicesDists.end(), IndexDist_Sorter());

				return nFound;
			}
		};

		//----------------------------------------------------------------------------------------------------------------------
		//			L2_Simple_Adaptor	
		//----------------------------------------------------------------------------------------------------------------------
		/** Squared Euclidean (L2) distance functor (suitable for low-dimensionality
		 * datasets, like 2D or 3D point clouds) Corresponding distance traits:
		 * nanoflann::metric_L2_Simple \tparam T Type of the elements (e.g. double,
		 * float, uint8_t) \tparam _DistanceType Type of distance variables (must be
		 * signed) (e.g. float, double, int64_t)
		 */
		template<typename T, typename DataSource, typename DT = T>
		struct L2_Simple_Adaptor {
			typedef T	ElementType;
			typedef DT	DistanceType;

			const DataSource& data_source;

			L2_Simple_Adaptor( const DataSource& ds ) : data_source( ds ){}

			inline DistanceType evalMetric(const T* a, const size_t b_idx, size_t size)const {
				DistanceType result = DistanceType();

				for (size_t idx = 0; idx < size; ++idx) {
					const DistanceType diff = a[idx] - data_source.kdtree_get_pt(b_idx, idx);
					result += diff * diff;
				}

				return result;
			}

			template<typename U, typename V>
			inline DistanceType accum_dist(const U a, const V b, const size_t)const {
				return (a - b) * (a - b);
			}
		};
	}// namespace KDTree
}// namespace dm