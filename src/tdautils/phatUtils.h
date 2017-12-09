#ifndef __PHATUTILS_H__
#define __PHATUTILS_H__

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

#include <utilities/timer.h>



template< typename Diagrams, typename PersistencePairs, typename Evaluator,
          typename Boundary >
inline void initDiagrams(
    Diagrams & diagrams, const PersistencePairs & pairs,
    const Evaluator & evaluator, const Boundary & boundary_matrix,
    const unsigned maxdimension) {

	diagrams.resize(maxdimension + 1);
	typename Diagrams::value_type::value_type dgmPoint(2);

	unsigned nPairs = pairs.get_num_pairs();

	// If persistence is not empty, manually add 0th homology for minimum to infinity
	if (nPairs > 0) {
		dgmPoint[0] = evaluator[0];
		dgmPoint[1] = std::numeric_limits< double >::infinity();
		diagrams[0].push_back(dgmPoint);
	}

  for (phat::index idx = 0; idx < nPairs; ++idx) {
    const phat::index b = pairs.get_pair(idx).first;
    const phat::index d = pairs.get_pair(idx).second;
    unsigned dim = (unsigned)boundary_matrix.get_dim(b);
		if (dim > maxdimension) {
			continue;
		}
		if (evaluator[b] < evaluator[d]) {
			dgmPoint[0] = evaluator[b];
			dgmPoint[1] = evaluator[d];
			diagrams[dim].push_back(dgmPoint);
		}
	}
}



template< typename Locations, typename PersistencePairs, typename Evaluator,
	        typename VectorList >
inline void initLocations(
    Locations & locations, const PersistencePairs & pairs,
    const Evaluator & evaluator, const VectorList & cmplx,
    const unsigned maxdimension) {

  unsigned verticesMax = 0;
  for (typename VectorList::const_iterator iVertices = cmplx.begin();
       iVertices != cmplx.end(); ++iVertices) {
    if (iVertices->size() == 1) {
      verticesMax = std::max(verticesMax, (unsigned)(*iVertices)[0]);
    }
  }

  // vertices range from 0 to verticesMax
  std::vector< double > verticesValues(
      verticesMax + 1, -std::numeric_limits< double >::infinity());

  unsigned iValue = 0;
  for (typename VectorList::const_iterator iVertices = cmplx.begin();
       iVertices != cmplx.end(); ++iVertices, ++iValue) {
    if (iVertices->size() == 1) {
      verticesValues[(*iVertices)[0]] = evaluator[iValue];
    }
  }

  typename Locations::value_type::value_type locPoint(2);
  locations.resize(maxdimension + 1);

  unsigned nPairs = pairs.get_num_pairs();

	// If persistence is not empty, manually trace back birth & death point
  if (nPairs > 0) {
    locPoint[0] = getLocation(cmplx[0], verticesValues);
    locPoint[1] = (unsigned)(
        std::max_element(verticesValues.begin(), verticesValues.end())
        - verticesValues.begin() + 1);
    locations[0].push_back(locPoint);
	}

  for (phat::index idx = 0; idx < nPairs; ++idx) {
    const phat::index b = pairs.get_pair(idx).first;
    const phat::index d = pairs.get_pair(idx).second;
    unsigned dim = cmplx[b].size() - 1;
    if (dim > maxdimension) {
			continue;
		}
    if (evaluator[b] < evaluator[d]) {
      locPoint[0] = getLocation(cmplx[b], verticesValues);
      locPoint[1] = getLocation(cmplx[d], verticesValues);
      locations[dim].push_back(locPoint);
		}
	}
}



// FiltrationDiag in PHAT
/** \brief Construct the persistence diagram from the filtration using library
 *         PHAT.
 *
 * @param[out] void             Void
 * @param[in]  cmplx            Simplicial complex of the input filtration
 * @param[in]  values           Filtration values of the input filtration
 * @param[in]  boundary_matrix  Boundary matrix of the input filtration
 * @param[in]  maxdimension     Max dimension of the homological features to be
 *                              computed
 * @param[in]  location         Are location of birth point, death point, and
 *                              representative cycles returned?
 * @param[in]  printProgress    Is progress printed?
 * @param[in]  persDgm          Memory space for the resulting persistence
 *                              diagram
 * @param[in]  persLoc          Memory space for the resulting birth points and
 *                              death points
 * @param[in]  persCycle        Memory space for the resulting representative
 *                              cycles
 */
template< typename VectorList, typename RealVector, typename Boundary >
void FiltrationDiagPhat(
  const VectorList                                      & cmplx,
  const RealVector                                      & values,
  Boundary                                              & boundary_matrix,
  const int                                               maxdimension,
  const bool                                              location,
  const bool                                              printProgress,
  std::vector< std::vector< std::vector< double > > >   & persDgm,
  std::vector< std::vector< std::vector< unsigned > > > & persLoc,
  std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {
	Timer persistence_timer;
	persistence_timer.start();

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
	initDiagrams(persDgm, pairs, values, boundary_matrix, maxdimension);

	// trace back birth & death simplex
	if (location) {
		initLocations(persLoc, pairs, values, cmplx, maxdimension);
	}

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}



# endif // __PHATUTILS_H__
