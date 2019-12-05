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
					};
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
			NodePtr					rootNode;
			size_t					leafMaxSize;
			size_t					size;				//!< Number of current points in the dataset
			size_t					sizeAtIndexBuild;	//!< Number of points in the dataset when the
														//!< index was built
			int						dim;				//!< Dimensionality of each data point
			BoundingBox				rootBox;			//!< The KD-tree used to find neighbours
			/**
			* Pooled memory allocator.
			*
			* Using a pooled memory allocator is more efficient
			* than allocating memory directly when there is a large
			* number small of memory allocations.
			*/
			PoolAllocator			pool;
		public:	// == FUNCTIONS ==
			/** Frees the previously-built index. Automatically called within
			  * buildIndex(). 
			  */
			void freeIndex(Derived &obj) {
				obj.pool.free_all();
				obj.rootNode = nullptr;
				obj.sizeAtIndexBuild = 0;
			}

			inline dataset_get(const Derived& obj, size_t idx, int dim)const {
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
				DistanceType splitVal = (box[cutfeat].low + box[cutfeat].high) / 2;
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

				if ((right - left) <= static_cast<IndexType>(obj.leafMaxSize)) {
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

					node->node_type.sub.divlow = left_bbox[cutfeat].high;
					node->node_type.sub.divHigh = right_box[cutfeat].low;
					
					for (int d = 0; d < DIM; d++) {
						bbox[d].low = std::min(left_bbox[d].low, right_box[d].low);
						bbox[d].high = std::max(left_bbox[d].high, right_box[d].high);
					}
				}

				return node;
			}
		};
	}
}