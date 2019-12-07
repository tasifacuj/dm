#pragma once

namespace dm {
	namespace dbscan {

	struct SearchParams {
		using Real = float;

		Real	Eps = 0.0f;	//!< search radius
		int		MinPts = 0;	//!< minimum number of points required to form a dense region

		SearchParams(Real epsilon, int minPts )
			: Eps( epsilon )
			, MinPts( minPts )
		{}
	};

		template< typename Dataset, typename DistanceSearch, typename IndexType>
		class DbScan {
		public: // == TYPES ==
			typedef typename DistanceSearch::ElementType ElementType;
		private:// == MEMBERS ==
			Dataset&				dataset_;
			const DistanceSearch&	searchFunc_;
		public:// == CTORs ==
			DbScan(DataSet& ds, const DistanceSearch& searchFunc)
				: dataset_( ds )
				, searchFunc_( searchFunc ){}

			void evaluate(const SearchParams& params) {
				int c = 0;
				for (IndexType idx = 0; idx < dataset_.size(); idx++) {
					if(  )
				}

			}
		};
	}
}
