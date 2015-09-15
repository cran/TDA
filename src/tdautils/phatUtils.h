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


template< typename Diagrams, typename PersistencePairs, typename Evaluator,
		typename SimplexMapInv >
inline void initDiagramsPhat(Diagrams& diagrams, const PersistencePairs& pairs,
		const Evaluator& evaluator, const SimplexMapInv simplex_map_inv,
		const unsigned maxdimension) {

	diagrams.resize(maxdimension + 1);
	typename Diagrams::value_type::value_type dgmPoint(2);

	// Manually add 0th homology for minimum to infinity
	dgmPoint[0] = evaluator(simplex_map_inv[0]);
	dgmPoint[1] = std::numeric_limits< double >::infinity();
	diagrams[0].push_back(dgmPoint);

	for (phat::index idx = 0; idx < pairs.get_num_pairs(); ++idx) {
		const typename SimplexMapInv::value_type& b =
				simplex_map_inv[pairs.get_pair(idx).first];
		const typename SimplexMapInv::value_type& d =
				simplex_map_inv[pairs.get_pair(idx).second];
		if ((unsigned)b.dimension() > maxdimension) {
			continue;
		}
		if (evaluator(b) < evaluator(d)) {
			dgmPoint[0] = evaluator(b);
			dgmPoint[1] = evaluator(d);
			diagrams[b.dimension()].push_back(dgmPoint);
		}
	}
}



template< typename Locations, typename PersistencePairs, typename Evaluator,
		typename SimplexMapInv, typename RealVector >
inline void initLocationsPhat(Locations& locations,
		const PersistencePairs& pairs, const Evaluator& evaluator, 
		const SimplexMapInv simplex_map_inv, const unsigned maxdimension,
		const RealVector& FUNvalues) {

	typename Locations::value_type::value_type locPoint(2);
	locations.resize(maxdimension + 1);

	// Manually trace back birth & death point
	locPoint[0] = getLocation(simplex_map_inv[0], FUNvalues);
	locPoint[1] = (unsigned)(std::max_element(
			FUNvalues.begin(), FUNvalues.end()) - FUNvalues.begin() + 1);
	locations[0].push_back(locPoint);

	for (phat::index idx = 0; idx < pairs.get_num_pairs(); ++idx) {
		const typename SimplexMapInv::value_type& b =
			simplex_map_inv[pairs.get_pair(idx).first];
		const typename SimplexMapInv::value_type& d =
			simplex_map_inv[pairs.get_pair(idx).second];
		if ((unsigned)b.dimension() > maxdimension) {
			continue;
		}
		if (evaluator(b) < evaluator(d)) {
			locPoint[0] = getLocation(b, FUNvalues);
			locPoint[1] = getLocation(d, FUNvalues);
			locations[b.dimension()].push_back(locPoint);
		}
	}
}



template< typename Fltr, typename Evaluator,
		typename RealVector >
void computePersistencePhat(Fltr f, const Evaluator& evaluator,
		const unsigned maxdimension, const RealVector& FUNvalues,
		const bool location, const bool printProgress,
		std::vector< std::vector< std::vector< double > > >& persDgm,
		std::vector< std::vector< std::vector< unsigned int > > >& persLoc) {

	Timer persistence_timer;
	persistence_timer.start();

	// convert from Dionysus to phat
	std::vector< typename Fltr::Simplex > simplex_map_inv;
	phat::boundary_matrix< phat::vector_vector > boundary_matrix;

	std::map< typename Fltr::Simplex, phat::index,
			typename Fltr::Simplex::VertexComparison > simplex_map;
	phat::index size_of_simplex_map = 0;
	simplex_map_inv.resize(f.size());
	boundary_matrix.set_num_cols(f.size());
	for (typename Fltr::Index it = f.begin(); it != f.end(); it++) {
		phat::column boundary_indices;
		const typename Fltr::Simplex& c = f.simplex(it);
		for(typename Fltr::Simplex::BoundaryIterator bit = c.boundary_begin();
				bit != c.boundary_end(); ++bit) {
			boundary_indices.push_back(simplex_map[*bit]);
		}
		std::sort(boundary_indices.begin(), boundary_indices.end());
		boundary_matrix.set_col(size_of_simplex_map, boundary_indices);
		phat::dimension dim_of_column =
				boundary_indices.size() == 0 ? 0 : boundary_indices.size() - 1;
		boundary_matrix.set_dim(size_of_simplex_map, dim_of_column);
		simplex_map_inv[size_of_simplex_map] = c;
		simplex_map.insert(typename std::map< typename Fltr::Simplex, phat::index >
				::value_type(c, size_of_simplex_map++));
	}

	// Compute persistent homology from sorted simplicial complex
	phat::persistence_pairs pairs;

	// choose an algorithm (choice affects performance)
	// and compute the persistence pair
	// (modifies boundary_matrix)
	phat::compute_persistence_pairs< phat::twist_reduction >(
			pairs, boundary_matrix);
	// sort the persistence pairs by birth index 
	pairs.sort();

	persistence_timer.stop();

	// Save persistent diagram
	initDiagramsPhat(persDgm, pairs, evaluator, simplex_map_inv, maxdimension);

	// trace back birth & death simplex
	if (location) {
		initLocationsPhat(persLoc, pairs, evaluator, simplex_map_inv,
				maxdimension, FUNvalues);
	}

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}
