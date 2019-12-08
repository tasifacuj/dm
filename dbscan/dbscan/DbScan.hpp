#pragma once
#include <algorithm>

namespace dm {
	namespace dbscan {

		struct ScanParams {
			using Real = float;

			Real	Eps = 0.0f;	//!< search radius
			int		MinPts = 0;	//!< minimum number of points required to form a dense region

			ScanParams(Real epsilon, int minPts)
				: Eps(epsilon)
				, MinPts(minPts)
			{}
		};

		template< typename DataSet, typename DistanceSearch, typename IndexType>
		class DbScan {
		public: // == TYPES ==
			typedef typename DistanceSearch::ElementType		ElementType;
			typedef typename  DistanceSearch::QueryResult		QueryResult;
			typedef typename DistanceSearch::QuerySearchParams	QuerySearchParams;
		private:// == MEMBERS ==
			DataSet&				dataset_;
			const DistanceSearch&	index_;
		public:// == CTORs ==
			DbScan(DataSet& ds, const DistanceSearch& dataIndex)
				: dataset_(ds)
				, index_(dataIndex) {}

			void evaluate(const ScanParams& scanParams) {
				int c = 0;																									/* Cluster counter */
				for (IndexType p_idx = 0; p_idx < dataset_.size(); p_idx++) {												/*for each point P in database DB*/
					if (dataset_.label(p_idx) != DataSet::Undefined)														/*if label(P) != undefined then continue*/
						continue;

					const ElementType search_radius = static_cast<ElementType>(scanParams.Eps);
					QueryResult  N_ret_matches;

					QuerySearchParams qparams;
					const auto& curPt = dataset_.getPt(p_idx);
					const ElementType query_pt[2] = { curPt.x, curPt.y };
					const size_t nMatches = index_.radiusSearch(&query_pt[0], search_radius, N_ret_matches, qparams);		/*Neighbors N = RangeQuery(DB, distFunc, P, eps)*/

					if (nMatches < scanParams.MinPts) {																		/*if |N| < minPts then */
						dataset_.setLabel(p_idx, DataSet::Noise);																/*label(P) = Noise*/
						continue;
					}

					c++;																									/* next cluster label */
					dataset_.setLabel(p_idx, c);																			/* label(P) = C */
					N_ret_matches.erase(p_idx);																				/* Seed set S = N \ {P}*/
					auto q_it = N_ret_matches.begin();

					for (auto q_it = N_ret_matches.begin(); q_it != N_ret_matches.end(); ) {								/*for each point Q in S*/
						if (dataset_.label(q_it->first) == DataSet::Noise)													/*if label(Q) = Noise then label(Q) = C*/
							dataset_.setLabel(q_it->first, c);

						if (dataset_.label(q_it->first) != DataSet::Undefined) {												/*if label(Q) != undefined then continue*/
							q_it++;
							continue;
						}

						dataset_.setLabel(q_it->first, c);																	/*label(Q) = C*/
						QuerySearchParams qparams2;
						const auto& curPt2 = dataset_.getPt(q_it->first);
						const ElementType query_pt2[2] = { curPt2.x, curPt2.y };
						QueryResult  N_ret_matches2;
						const size_t nMatches2 = index_.radiusSearch(&query_pt2[0], search_radius, N_ret_matches2, qparams2);	/*Neighbors N = RangeQuery(DB, distFunc, Q, eps)*/

						if (nMatches2 > scanParams.MinPts) {																	/*if |N| >= minPts then S = S | N */
							N_ret_matches.insert(N_ret_matches2.begin(), N_ret_matches2.end());
							q_it = N_ret_matches.begin();
						}else{
							q_it++;
						}
					}

				}
			}
		};
	} // namespace dbscan
}// namespace dm
