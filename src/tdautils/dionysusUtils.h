#ifndef __DIONYSUSUTILS_H__
#define __DIONYSUSUTILS_H__

#include <topology/simplex.h>
#include <string>

#include <utilities/timer.h>



typedef std::vector< double > Point;
typedef std::vector< Point > PointContainer;



template< typename Simplex >
class ComparisonDataDimension :
  public std::binary_function<const Simplex &, const Simplex &, bool> {

public:
  bool operator()(const Simplex & a, const Simplex & b) const {
    if (a.data() == b.data())
      return a.dimension() < b.dimension();
    else
      return a.data() < b.data();
  }
};



/**
 * Class: EvaluatePushBack<Container>
 *
 * Push back the simplex and the evaluated value
 */
template< typename Container, typename Evaluator >
class EvaluatePushBack {

public:
  EvaluatePushBack(Container & argContainer, const Evaluator & argEvaluator) :
    container(argContainer), evaluator(argEvaluator) {}

  void operator()(const typename Container::value_type & argSmp) const {
    typename Container::value_type smp(argSmp.vertices(), evaluator(argSmp));
    container.push_back(smp);
  }

private:
  Container & container;
  const Evaluator & evaluator;
};



template< typename VertexList, typename Evaluator >
unsigned getLocation(const VertexList & vertices, const Evaluator & evaluator) {
  typename VertexList::const_iterator vertexItr;
  unsigned vertex = *(vertices.begin());
	for (vertexItr = vertices.begin(); vertexItr != vertices.end(); ++vertexItr) {
		if (evaluator[*vertexItr] > evaluator[vertex]) {
			vertex = *vertexItr;
		}
	}
	return vertex + 1;
}



template< typename Diagrams, typename Iterator, typename Evaluator,
          typename SimplexMap >
inline void initDiagrams(
    Diagrams & diagrams, const Iterator & bg, const Iterator & end,
    const Evaluator & evaluator, const SimplexMap & m,
		const unsigned maxdimension) {

	diagrams.resize(maxdimension + 1);
	typename Diagrams::value_type::value_type dgmPoint(2);
	for (Iterator cur = bg; cur != end; ++cur) {
		// positive simplices corresponds to
		// negative simplices having non-empty cycles
		if (cur->sign()) {
			// first consider that cycle is paired
			if (!cur->unpaired()) {
				// the cycle that was born at cur is killed 
				// when we added death (another simplex)
				const typename SimplexMap::key_type & death = cur->pair;

				const typename SimplexMap::value_type & b = m[cur];
				const typename SimplexMap::value_type & d = m[death];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				if (evaluator(b) < evaluator(d)) {
					dgmPoint[0] = evaluator(b);
					dgmPoint[1] = evaluator(d);
					diagrams[b.dimension()].push_back(dgmPoint);
				}
			}
			else {    // cycles can be unpaired
				const typename SimplexMap::value_type& b = m[cur];
				dgmPoint[0] = evaluator(b);
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				dgmPoint[1] = std::numeric_limits< double >::infinity();
				diagrams[b.dimension()].push_back(dgmPoint);
			}
		}
	}
}



template< typename Simplex, typename Locations, typename Cycles,
          typename Persistence, typename Evaluator, typename SimplexMap,
          typename Filtration >
inline void initLocations(
    Locations & locations, Cycles & cycles, const Persistence & p,
    const Evaluator & evaluator, const SimplexMap & m,
    const unsigned maxdimension, const Filtration & filtration) {

	unsigned verticesMax = 0;
	for (typename Filtration::Index iFltr = filtration.begin();
       iFltr != filtration.end(); ++iFltr) {
		const typename Filtration::Simplex & c = filtration.simplex(iFltr);
		if (c.dimension() == 0) {
			verticesMax = std::max(verticesMax, *(c.vertices().begin()));
		}
	}

  // vertices range from 0 to verticesMax
  std::vector< double > verticesValues(
      verticesMax + 1, -std::numeric_limits< double >::infinity());

  for (typename Filtration::Index iFltr = filtration.begin();
       iFltr != filtration.end(); ++iFltr) {
		const typename Filtration::Simplex & c = filtration.simplex(iFltr);
		if(c.dimension() == 0) {
			verticesValues[*(c.vertices().begin())] = c.data();
		}
	}

	locations.resize(maxdimension + 1);
	cycles.resize(maxdimension + 1);
	typename Locations::value_type::value_type persLocPoint(2);
	typename Cycles::value_type::value_type persBdy;
	typename Cycles::value_type::value_type::value_type persSimplex;
	for (typename Persistence::iterator cur = p.begin(); cur != p.end(); ++cur) {
		// positive simplices corresponds to
		// negative simplices having non-empty cycles
		if (cur->sign()) {
			// first consider that cycle is paired
			if (!cur->unpaired()) {
				// the cycle that was born at cur is killed 
				// when we added death (another simplex)
				const typename SimplexMap::key_type& death = cur->pair;

				//const typename SimplexMap::value_type& b = m[cur];
				//const typename SimplexMap::value_type& d = m[death];
        const typename Filtration::Simplex & b = m[cur];
        const typename Filtration::Simplex & d = m[death];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				if (evaluator(b) < evaluator(d)) {
					persLocPoint[0] = getLocation(b.vertices(), verticesValues);
					persLocPoint[1] = getLocation(d.vertices(), verticesValues);
					locations[b.dimension()].push_back(persLocPoint);

					// Iterate over the cycle
					persBdy.clear();
					const typename Persistence::Cycle& cycle = death->cycle;
					for (typename Persistence::Cycle::const_iterator
						si = cycle.begin(); si != cycle.end(); ++si) {
						persSimplex.clear();
						const typename Simplex::VertexContainer&
							vertices = m[*si].vertices();    // std::vector<Vertex> where Vertex = Distances::IndexType
						typename Simplex::VertexContainer::const_iterator vtxItr;
						for (vtxItr = vertices.begin(); vtxItr != vertices.end();
							++vtxItr) {
							persSimplex.push_back(*vtxItr + 1);
						}
						persBdy.push_back(persSimplex);
					}
					cycles[b.dimension()].push_back(persBdy);
				}
			}
			else {    // cycles can be unpaired
				const typename SimplexMap::value_type& b = m[cur];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				persLocPoint[0] = getLocation(b.vertices(), verticesValues);
				persLocPoint[1] = (unsigned)(
            std::max_element(verticesValues.begin(), verticesValues.end())
            - verticesValues.begin() + 1);
				locations[b.dimension()].push_back(persLocPoint);

				// Iterate over the cycle
				persBdy.clear();
				cycles[b.dimension()].push_back(persBdy);
			}
		}
	}
}



// FiltrationDiag in Dionysus
/** \brief Construct the persistence diagram from the filtration using library
 *         Dionysus.
 *
 * @param[out] void           Void
 * @param[in]  filtration     The input filtration
 * @param[in]  maxdimension   Max dimension of the homological features to be
 *                            computed
 * @param[in]  location       Are location of birth point, death point, and
 *                            representative cycles returned?
 * @param[in]  printProgress  Is progress printed?
 * @param[in]  persDgm        Memory space for the resulting persistence
 *                            diagram
 * @param[in]  persLoc        Memory space for the resulting birth points and
 *                            death points
 * @param[in]  persCycle      Memory space for the resulting representative
 *                            cycles
 * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
 *                            diagram. Diagram must point to enough memory space for
 *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
 *                            write nothing after.
 */
template< typename Persistence, typename Filtration >
void FiltrationDiagDionysus(
    const Filtration                                      & filtration,
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
	Persistence p(filtration); // initialize persistence
	if (!location) {
		p.pair_simplices(printProgress); // pair simplices
	}
	else {
		if (printProgress) {
			p.pair_simplices(p.begin(), p.end(), true,
				typename Persistence::PairVisitor(p.size()));
		}
		else {
			p.pair_simplices(p.begin(), p.end(), true,
				typename Persistence::PairVisitorNoProgress());
		}
	}

	persistence_timer.stop();

	// Save persistent diagram
  typename Persistence::template SimplexMap< Filtration >
      m = p.make_simplex_map(filtration);
  initDiagrams(persDgm, p.begin(), p.end(),
      typename Filtration::Simplex::DataEvaluator(), m, maxdimension);

	// TODO: why doesn't this work? rLog(rlmain, "testing");  
	//Persistence::SimplexMap< Fltr > m = p.make_simplex_map(f);
	//std::map<Dimension, PersistenceDiagram<> > dgms;
	//init_diagrams(dgms, p.begin(), p.end(),
	//	evaluate_through_map(m, Smplx::DataEvaluator()),
	//	evaluate_through_map(m, Smplx::DimensionExtractor()));

//#if 1
//	// Output cycles
//	PersistenceR::SimplexMap<FltrR>   m = p.make_simplex_map(f);
//	for (PersistenceR::iterator cur = p.begin(); cur != p.end(); ++cur)
//	{
//		//			const PersistenceR::Cycle& cycle = cur->cycle;
//
//		if (!cur->sign())        // only negative simplices have non-empty cycles
//		{
//			PersistenceR::OrderIndex birth = cur->pair;      // the cycle that cur killed was born when we added birth (another simplex)
//
//			const SmplxR& b = m[birth];
//			const SmplxR& d = m[cur];
//
//			// if (b.dimension() != 1) continue;
//			// std::cout << "Pair: (" << size(b) << ", " << size(d) << ")" << std::endl;
//			if (b.dimension() > maxdimension) continue;
//			diagram_out << b.dimension() << " " << size(b) << " " << size(d) << std::endl;
//		}
//		else if (cur->unpaired())    // positive could be unpaired
//		{
//			const SmplxR& b = m[cur];
//			// if (b.dimension() != 1) continue;
//
//			// std::cout << "Unpaired birth: " << size(b) << std::endl;
//			// cycle = cur->chain;      // TODO
//			if (b.dimension() > maxdimension) continue;
//			diagram_out << b.dimension() << " " << size(b) << " inf" << std::endl;
//		}
//
//		// Iterate over the cycle
//		// for (PersistenceR::Cycle::const_iterator si =  cycle.begin();
//		//                                                          si != cycle.end();     ++si)
//		// {
//		//     const SmplxR& s = m[*si];
//		//     //std::cout << s.dimension() << std::endl;
//		//     const SmplxR::VertexContainer& vertices = s.vertices();          // std::vector<Vertex> where Vertex = Distances::IndexType
//		//     AssertMsg(vertices.size() == s.dimension() + 1, "dimension of a simplex is one less than the number of its vertices");
//		//     std::cout << vertices[0] << " " << vertices[1] << std::endl;
//		// }
//	}
//#endif

	// trace back birth & death simplex
	if (location) {
    initLocations< typename Filtration::Simplex >(
        persLoc, persCycle, p, typename Filtration::Simplex::DataEvaluator(),
        m, maxdimension, filtration);
	}

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}



// RipsFiltration in Dionysus
/** \brief Construct the Rips filtration constructed on the input set of points
 *         using library GUDHI.
 *
 * @param[out] Filtration     A filtration
 * @param[in]  X              Either an nxd matrix of coordinates,
 *                            or an nxn matrix of distances of points
 * @param[in]  is_row_names   Wehther row names are included in the input X
 * @param[in]  maxdimension   Max dimension of the homological features to be
 *                            computed.
 * @param[in]  maxscale       Threshold for the Rips complex
 * @param[in]  printProgress  Is progress printed?
 * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
 *                            diagram. Diagram must point to enough memory space for
 *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
 *                            write nothing after.
 */
template< typename Distances, typename Generator, typename Filtration,
          typename RealMatrix, typename Print >
inline Filtration RipsFiltrationDionysus(
    const RealMatrix & X,
    const unsigned     nSample, 
    const unsigned     nDim,
    const bool         is_row_names,
    const int          maxdimension,
    const double       maxscale,
    const bool         printProgress,
    const Print      & print
) {

  PointContainer points = TdaToStl< PointContainer >(X, nSample, nDim,
      is_row_names);
  //read_points(infilename, points);
  //read_points2(infilename, points);

  Distances distances(points);
  Generator rips(distances);
  typename Generator::Evaluator size(distances);
  Filtration filtration;
  EvaluatePushBack< Filtration, typename Generator::Evaluator > functor(
      filtration, size);

  // Generate maxdimension skeleton of the Rips complex
  rips.generate(maxdimension + 1, maxscale, functor);

  if (printProgress) {
    print("# Generated complex of size: %d \n", filtration.size());
  }

  // Sort the simplices with respect to comparison criteria
  // e.g. distance or function values
  filtration.sort(ComparisonDataDimension< typename Filtration::Simplex >());

  return filtration;
}



# endif // __DIONYSUSUTILS_H__
