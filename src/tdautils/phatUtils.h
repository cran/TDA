// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

#include <phat/boundary_matrix.h>

#include <limits>
#include <algorithm>

#include <tdautils/dionysusUtils.h>


template<typename Flt>
void computePersistentPairsPhat(Flt f, int maxDimension, const double * const FUNvaluesInput, const unsigned int & gridNumProd, bool isLocation, std::vector< std::vector< std::vector< double > > > &persDgm, std::vector< std::vector< std::vector< unsigned int > > > &persLoc)
{

	// If phat is used, convert from Dionysus to phat
	std::vector<Smplx> simplex_map_inv;
	phat::boundary_matrix< phat::vector_vector > boundary_matrix;

	std::map<Smplx,phat::index,Smplx::VertexComparison> simplex_map;
	phat::index size_of_simplex_map=0;
	simplex_map_inv.resize( f.size() );
	boundary_matrix.set_num_cols( f.size() );
	for(Fltr::Index it=f.begin();it!=f.end();it++)
	{
		phat::column boundary_indices;
		const Smplx& c = f.simplex(it);
		for(Smplx::BoundaryIterator bit = c.boundary_begin(); bit != c.boundary_end(); bit++)
		{
			boundary_indices.push_back( simplex_map[*bit] );
		}
		std::sort(boundary_indices.begin(),boundary_indices.end());
		boundary_matrix.set_col( size_of_simplex_map, boundary_indices );
		phat::dimension dim_of_column = boundary_indices.size()==0 ? 0 : boundary_indices.size()-1;
		boundary_matrix.set_dim( size_of_simplex_map, dim_of_column );
		simplex_map_inv[size_of_simplex_map] = c;
		simplex_map.insert(std::map<Smplx,phat::index>::value_type(c, size_of_simplex_map++));
	}

	// Compute persistent homology from sorted simplicial complex
	phat::persistence_pairs pairs;

	// choose an algorithm (choice affects performance) and compute the persistence pair
	// (modifies boundary_matrix)
	phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );
	// sort the persistence pairs by birth index 
	pairs.sort();


	// Save persistent diagram & trace back birth & death simplex
	std::vector< double > persDgmPoint(2);
	std::vector< unsigned int > persLocPoint(2);
	unsigned int persDim;

	// Manually add 0th homology for minimum to infinity
	persDgmPoint[0] = simplex_map_inv.at(0).data();
	persDgmPoint[1] = std::numeric_limits< double >::infinity();
	persDgm[0].push_back(persDgmPoint );

	// Manually trace back birth & death point
	if (isLocation)
	{
		persLocPoint[0] = getLocation(simplex_map_inv.at(0), FUNvaluesInput);
		persLocPoint[1] = (unsigned int)(std::max_element(FUNvaluesInput, FUNvaluesInput+gridNumProd)-FUNvaluesInput+1);
		persLoc[ 0 ].push_back( persLocPoint );
	}

	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
		persDim = simplex_map_inv.at(pairs.get_pair( idx ).first).dimension();
		persDgmPoint[0] = simplex_map_inv.at(pairs.get_pair( idx ).first).data();
		persDgmPoint[1] = simplex_map_inv.at(pairs.get_pair( idx ).second).data();
		if (persDgmPoint[0] < persDgmPoint[1] && persDim <= maxDimension )
		{
			persDgm[ persDim ].push_back( persDgmPoint );

			// trace back birth & death point
			if (isLocation)
			{
				persLocPoint[0] = getLocation(simplex_map_inv.at(pairs.get_pair( idx ).first), FUNvaluesInput);
				persLocPoint[1] = getLocation(simplex_map_inv.at(pairs.get_pair( idx ).second), FUNvaluesInput);
				persLoc[ persDim ].push_back( persLocPoint );
			}
		}
	}
}
