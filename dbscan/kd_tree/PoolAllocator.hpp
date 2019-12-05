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
#pragma once
#include <cstdlib>

namespace dm {
	namespace KdTree {

		/**
		 * Pooled storage allocator
		 *
		 * The following routines allow for the efficient allocation of storage in
		 * small chunks from a specified pool.  Rather than allowing each structure
		 * to be freed individually, an entire pool of storage is freed at once.
		 * This method has two advantages over just using malloc() and free().  First,
		 * it is far more efficient for allocating small objects, as there is
		 * no overhead for remembering all the information needed to free each
		 * object or consolidating fragmented memory.  Second, the decision about
		 * how long to keep an object is made at the time of allocation, and there
		 * is no need to track down all the objects to free them.
		 *
		 * We maintain memory alignment to word boundaries by requiring that all allocations be in multiples of the machine wordsize.
	     * Size of machine word in bytes.  Must be power of 2.
	     * Minimum number of bytes requested at a time from	the system.  Must be
	     * multiple of WORDSIZE. 
		 */

		class PoolAllocator {
		public: // == CONSTANTS ==
			static constexpr size_t BlockSize = 8192;
			static constexpr size_t WordSize = 16;
		private:// == MEMBERS ==
			size_t remaining_;	//!< Number of bytes left in current block of storage.
			void* base_;		//!< Pointer to base of current block of storage.
			void* loc_;			//!< Current location in block to next allocate memory.
		public:
			size_t usedMemory_;
			size_t wastedMemory_;
		public: // == CTORs ==
			PoolAllocator() { internalInit(); }
			~PoolAllocator() { freeAll(); }
		public:
			void freeAll() {
				while (base_ != nullptr) {
					void* prev = *(static_cast<void**>(base_));
					::free(base_);
					base_ = prev;
				}

				internalInit();
			}

			/**
			 * Allocates (using this pool) a generic type T.
			 *
			 * Params:
			 *     count = number of instances to allocate.
			 * Returns: pointer (of type T*) to memory buffer
			 */
			template<typename T>
			inline T* allocate(const size_t count = 1) {
				T* mem = static_cast<T*>(this->mallocPrivate(sizeof(T) * count));
				return mem;
			}
		private:
			void internalInit() {
				remaining_ = 0;
				base_ = nullptr;
				usedMemory_ = 0;
				wastedMemory_ = 0;
			}

			/**
			 * Returns a pointer to a piece of new memory of the given size in bytes
			 * allocated from the pool.
			 */
			void* mallocPrivate(const size_t req_size) {
				/* Round size up to a multiple of wordsize.  The following expression
				   only works for WORDSIZE that is a power of 2, by masking last bits of
				   incremented size to zero.
				*/
				const size_t size = (req_size + (WordSize - 1)) & ~(WordSize - 1);
				/* Check whether a new block must be allocated.  Note that the first word
					of a block is reserved for a pointer to the previous block.
				 */

				if (size > remaining_) {
					wastedMemory_ += remaining_;
					// allocate new storage
					const size_t blocksize = (size + sizeof(void*) + (WordSize - 1) > BlockSize) ? size + sizeof(void*) + (WordSize - 1) : BlockSize;
					void* m = ::malloc(blocksize);

					if (!m)
						return nullptr;

					/* Fill first word of new block with pointer to previous block. */
					static_cast<void**>(m)[0] = base_;
					base_ = m;
					size_t shift = 0;
					remaining_ = blocksize - sizeof(void*) - shift;
					loc_ = (static_cast<char*>(m) + sizeof(void*) + shift);
				}

				void* rloc = loc_;
				loc_ = static_cast<char*>(loc_) + size;
				remaining_ -= size;
				usedMemory_ += size;
				return rloc;
			}
		};
	}// namespace KdTree
}// namespace dm
