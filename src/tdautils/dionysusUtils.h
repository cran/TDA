#ifndef __DIONYSUSUTILS_H__
#define __DIONYSUSUTILS_H__

#include <topology/simplex.h>
#include <string>
#include <sstream>
#include <cstdlib>



template< typename Simplex >
std::vector< unsigned > getVertices(const Simplex& smp) {
	std::stringstream sstr;
	std::vector<unsigned int> vertices(smp.dimension() + 1);
	sstr << smp;
	std::string vtxStr;
	std::getline(sstr,vtxStr,'<');
	unsigned int vtxIdx;
	for (vtxIdx = 0; vtxIdx < (unsigned)smp.dimension(); ++vtxIdx)
	{
		std::getline(sstr, vtxStr, ',');
		vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	}
	std::getline(sstr,vtxStr,'>');
	vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	return vertices;
}



template< typename Simplex, typename Evaluator >
unsigned getLocation(const Simplex& smp, const Evaluator& evaluator) {
	std::vector< unsigned > vertices;
	std::vector< unsigned >::const_iterator vertexItr;
	unsigned vertex;
	vertices = getVertices(smp);
	vertex = *(vertices.begin());
	for (vertexItr = vertices.begin(); vertexItr != vertices.end(); ++vertexItr)
	{
		if (evaluator[*vertexItr] > evaluator[vertex])
		{
			vertex = *vertexItr;
		}
	}
	return vertex + 1;
}



template< typename Diagrams, typename Iterator, typename Evaluator,
		typename SimplexMap >
inline void initDiagrams(Diagrams& diagrams, const Iterator& bg,
		const Iterator& end, const Evaluator& evaluator, const SimplexMap& m,
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
				const typename SimplexMap::key_type& death = cur->pair;

				const typename SimplexMap::value_type& b = m[cur];
				const typename SimplexMap::value_type& d = m[death];
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
		typename RealVector >
inline void initLocations(Locations& locations, Cycles& cycles,
		const Persistence& p, const Evaluator& evaluator, const SimplexMap& m,
		const unsigned maxdimension, const RealVector& FUNvalues) {

	locations.resize(maxdimension + 1);
	cycles.resize(maxdimension + 1);
	typename Locations::value_type::value_type persLocPoint(2);
	typename Cycles::value_type::value_type persCyclePoint;
	for (typename Persistence::iterator cur = p.begin(); cur != p.end(); ++cur) {
		// positive simplices corresponds to
		// negative simplices having non-empty cycles
		if (cur->sign()) {
			// first consider that cycle is paired
			if (!cur->unpaired()) {
				// the cycle that was born at cur is killed 
				// when we added death (another simplex)
				const typename SimplexMap::key_type& death = cur->pair;

				const typename SimplexMap::value_type& b = m[cur];
				const typename SimplexMap::value_type& d = m[death];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				if (evaluator(b) < evaluator(d)) {
					persLocPoint[0] = getLocation(b, FUNvalues);
					persLocPoint[1] = getLocation(d, FUNvalues);
					locations[b.dimension()].push_back(persLocPoint);

					// Iterate over the cycle
					persCyclePoint.clear();
					const typename Persistence::Cycle& cycle = death->cycle;
					for (typename Persistence::Cycle::const_iterator
						si = cycle.begin(); si != cycle.end(); ++si) {
						const typename Simplex::VertexContainer&
								vertices = m[*si].vertices();    // std::vector<Vertex> where Vertex = Distances::IndexType
						typename Simplex::VertexContainer::const_iterator vtxItr;
						for (vtxItr = vertices.begin(); vtxItr != vertices.end();
						++vtxItr) {
							persCyclePoint.insert(*vtxItr + 1);
						}
					}
					cycles[b.dimension()].push_back(persCyclePoint);
				}
			}
			else {    // cycles can be unpaired
				const typename SimplexMap::value_type& b = m[cur];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				persLocPoint[0] = getLocation(b, FUNvalues);
				persLocPoint[1] = (unsigned)(std::max_element(
					FUNvalues.begin(), FUNvalues.end()) - FUNvalues.begin() + 1);
				locations[b.dimension()].push_back(persLocPoint);

				// Iterate over the cycle
				persCyclePoint.clear();
				cycles[b.dimension()].push_back(persCyclePoint);
			}
		}
	}
}



template< typename Persistence, typename Fltr, typename Evaluator,
		typename RealVector >
void computePersistenceDionysus(Fltr f, const Evaluator& evaluator,
		const unsigned maxdimension, const RealVector& FUNvalues,
		const bool location, const bool printProgress,
		std::vector< std::vector< std::vector< double > > > &persDgm,
		std::vector< std::vector< std::vector< unsigned int > > > &persLoc,
		std::vector< std::vector< std::set< unsigned int > > >& persCycle) {

	Timer persistence_timer;
	persistence_timer.start();

	// Compute persistent homology from sorted simplicial complex
	Persistence p(f); // initialize persistence
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
	typename Persistence::template SimplexMap< Fltr >
			m = p.make_simplex_map(f);
	initDiagrams(persDgm, p.begin(), p.end(), evaluator, m, maxdimension);

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
		initLocations< typename Fltr::Simplex >
				(persLoc, persCycle, p, evaluator, m, maxdimension, FUNvalues);
	}

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}



# endif // __DIONYSUSUTILS_H__
