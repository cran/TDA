// for R
#include <R.h>
#include <R_ext/Print.h>

// for Rcpp
#include <Rcpp.h>

// for Rips
#include <tdautils/ripsL2.h>
#include <tdautils/ripsArbit.h>

// for grid
#include <tdautils/gridUtils.h>

// for kernel density
#include <tdautils/kernelUtils.h>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

//for CGAL
#include <tdautils/cgalUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>



// GridDiag by Brittany T. Fasy
// modified by Jisu Kim for
// arbitrary dimension & using memory as an input & setting maximum dimension.
/** \brief Interface for R code, construct the persistence diagram
  * of sublevel/superlevel sets of a function evaluated over a grid of points
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  FUNvalues      A vector of length m1*...*md of function values over grid
  * @param[in]  gridDim        A vector (m1, ..., md),
  *                            where mi is number of grid points in ith dimension
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  decomposition  Either "5tetrahedra" or "barycenter"
  * @param[in]  library        Either "Dionysus" or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  */
 // [[Rcpp::export]]
Rcpp::List
GridDiag(const Rcpp::NumericVector & FUNvalues
       , const Rcpp::IntegerVector & gridDim
       , const int                   maxdimension
       , const std::string         & decomposition
       , const std::string         & library
       , const bool                  location
       , const bool                  printProgress
	) {
#ifdef LOGGING
	//rlog::RLogInit(argc, argv);

	stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
#endif

	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	Fltr f;

	// Generate simplicial complex from function values and grid
	if (decomposition[0] == '5') {
		simplicesFromGrid(f, FUNvalues, gridDim, maxdimension + 1);
	}
	if (decomposition[0] == 'b') {
		simplicesFromGridBarycenter(f, FUNvalues, gridDim, maxdimension + 1);
	}
	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", f.size());
	}

	// Sort the simplices with respect to function values
	f.sort(Smplx::DataComparison());

	// Compute the persistence diagram of the complex
	if (library[0] == 'D') {
		computePersistenceDionysus< Persistence >(f, Smplx::DataEvaluator(),
				maxdimension, FUNvalues, location, printProgress,
				persDgm, persLoc, persCycle);
	}
	if (library[0] == 'P') {
		 computePersistencePhat(f, Smplx::DataEvaluator(), maxdimension,
				FUNvalues, location, printProgress, persDgm, persLoc);
	}

	// Output persistent diagram
	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
		concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
		StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// RipsDiag for L2 distance in Dionysus or PHAT
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nxd matrix of coordinates,
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  library        Either "Dionysus" or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  */
template< typename RealMatrix >
inline Rcpp::List
RipsDiagL2DionysusPhat(const RealMatrix  & X 
                     , const int           maxdimension
                     , const double        maxscale
                     , const std::string & library
                     , const bool          location
	                 , const bool          printProgress
	) {
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;
	
	PointContainer points = RcppToStl< PointContainer >(X);
	//read_points(infilename, points);

	PairDistances           distances(points);
	Generator               rips(distances);
	Generator::Evaluator    size(distances);
	FltrR                   f;

	// Generate 2-skeleton of the Rips complex for epsilon = 50
	rips.generate(maxdimension + 1, maxscale, make_push_back_functor(f));

	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", f.size());
	}

	// Sort the simplices with respect to distance 
	f.sort(Generator::Comparison(distances));

	// Compute the persistence diagram of the complex
	if (library[0] == 'D') {
		computePersistenceDionysus< PersistenceR >(f, size, maxdimension,
				Rcpp::NumericVector(X.nrow()), location, printProgress,
				persDgm, persLoc, persCycle);
	}
	if (library[0] == 'P') {
		computePersistencePhat(f, size, maxdimension,
				Rcpp::NumericVector(X.nrow()), location, printProgress,
				persDgm, persLoc);
	}

	// Output persistent diagram
	return Rcpp::List::create(
			concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
			concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
			StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// RipsDiag for arbitrary distance in Dionysus or PHAT
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  library        Either "Dionysus" or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  */
template< typename RealMatrix >
inline Rcpp::List
RipsDiagArbitDionysusPhat(const RealMatrix  & X
                        , const int           maxdimension
                        , const double        maxscale
                        , const std::string & library
                        , const bool          location
                        , const bool          printProgress
	) {
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	PointContainer points = RcppToStl< PointContainer >(X, true);
	//read_points2(infilename, points);

	PairDistancesA           distances(points);
	GeneratorA               rips(distances);
	GeneratorA::Evaluator    size(distances);
	FltrRA                   f;

	// Generate 2-skeleton of the Rips complex for epsilon = 50
	rips.generate(maxdimension + 1, maxscale, make_push_back_functor(f));

	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", f.size());
	}

	// Sort the simplices with respect to distance 
	f.sort(GeneratorA::Comparison(distances));

	// Compute the persistence diagram of the complex
	if (library[0] == 'D') {
		computePersistenceDionysus< PersistenceR >(f, size, maxdimension,
				Rcpp::NumericVector(X.nrow()), location, printProgress,
				persDgm, persLoc, persCycle);
	}
	if (library[0] == 'P') {
		computePersistencePhat(f, size, maxdimension,
				Rcpp::NumericVector(X.nrow()), location, printProgress,
				persDgm, persLoc);
	}

	// Output persistent diagram
	return Rcpp::List::create(
			concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
			concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
			StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// [[Rcpp::export]]
double
Bottleneck(const Rcpp::NumericMatrix & Diag1
         , const Rcpp::NumericMatrix & Diag2
	) {
	return bottleneck_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2));
}



// [[Rcpp::export]]
double
Wasserstein(const Rcpp::NumericMatrix & Diag1
          , const Rcpp::NumericMatrix & Diag2
          , const int                   p
	) {
	return wasserstein_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2), p);
}



// KDE function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
Kde(const Rcpp::NumericMatrix & X
  , const Rcpp::NumericMatrix & Grid
  , const double                h
  , const Rcpp::NumericVector & weight
  , const bool                  printProgress
	) {
	const double pi = 3.141592653589793;
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	const double den = pow(h, (int)dimension) * pow(2 * pi, dimension / 2.0);
	Rcpp::NumericVector kdeValue;
	int counter = 0, percentageFloor = 0;
	int totalCount = gridNum;

	if (printProgress) {
		printProgressFrame(Rprintf);
	}

	if (dimension <= 1) {
		kdeValue = computeKernel< Rcpp::NumericVector >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}
	else {
		kdeValue = computeGaussOuter< Rcpp::NumericVector >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);

	}

	for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
		kdeValue[gridIdx] /= den;
	}

	if (printProgress) {
		Rprintf("\n");
	}

	return kdeValue;
}



// kernel Dist function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
KdeDist(const Rcpp::NumericMatrix & X
      , const Rcpp::NumericMatrix & Grid
      , const double                h
      , const Rcpp::NumericVector & weight
      , const bool printProgress
	) {
	const unsigned sampleNum = X.nrow();
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	// first = sum K_h(X_i, X_j), second = K_h(x, x), third = sum K_h(x, X_i)
	std::vector< double > firstValue;
	const double second = 1.0;
	std::vector< double > thirdValue;
	double firstmean;
	Rcpp::NumericVector kdeDistValue(gridNum);
	int counter = 0, percentageFloor = 0;
	int totalCount = sampleNum + gridNum;

	if (printProgress) {
		printProgressFrame(Rprintf);
	}

	firstValue = computeKernel< std::vector< double > >(
			X, X, h, weight, printProgress, Rprintf, counter, totalCount,
			percentageFloor);

	if (dimension <= 1) {
		thirdValue = computeKernel< std::vector< double > >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}
	else {
		thirdValue = computeGaussOuter< std::vector< double > >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}

	if (weight.size() == 1) {
		firstmean = std::accumulate(firstValue.begin(), firstValue.end(), 0.0) / sampleNum;
	}
	else {
		firstmean = std::inner_product(
				firstValue.begin(), firstValue.end(), weight.begin(), 0.0) / 
				std::accumulate(weight.begin(), weight.end(), 0.0);
	}

	for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
		kdeDistValue[gridIdx] = std::sqrt(firstmean + second - 2 * thirdValue[gridIdx]);
	}

	if (printProgress) {
		Rprintf("\n");
	}

	return kdeDistValue;
}



// distance to measure function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
Dtm(const Rcpp::NumericMatrix & knnDistance
  , const double                weightBound
  , const double                r
	) {
	const unsigned gridNum = knnDistance.nrow();
  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
	Rcpp::NumericVector dtmValue(gridNum, 0.0);
  unsigned weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += distanceTemp * distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += std::pow(distanceTemp, r);
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}



// distance to measure function on a Grid, with weight
// [[Rcpp::export]]
Rcpp::NumericVector
DtmWeight(const Rcpp::NumericMatrix & knnDistance
        , const double                weightBound
        , const double                r
        , const Rcpp::NumericMatrix & knnIndex
        , const Rcpp::NumericVector & weight
  ) {
  const unsigned gridNum = knnDistance.nrow();
  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
  Rcpp::NumericVector dtmValue(gridNum, 0.0);
  double weightTemp, weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += distanceTemp * distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
          (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * (weightBound - weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += std::pow(distanceTemp, r) * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
          (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}



// RipsDiag in GUDHI
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  printProgress  Is progress printed?
  * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
  *                            diagram. Diagram must point to enough memory space for
  *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
  *                            write nothing after.
  */
inline Rcpp::List
RipsDiagGUDHI(const Rcpp::NumericMatrix & X          //points to some memory space
            , const int                   maxdimension
            , const double                maxscale
            , const bool                  printProgress
	) {
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	// Turn the input points into a range of points
	typedef std::vector< double > Point_t;
	std::vector< Point_t > point_set =
			RcppToStl< std::vector< Point_t > >(X);


	// Compute the proximity graph of the points
	Graph_t prox_graph = compute_proximity_graph(point_set, maxscale
		, euclidean_distance<Point_t>);

	// Construct the Rips complex in a Simplex Tree
	Gudhi::Simplex_tree<> st;
	st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
	st.expansion(maxdimension + 1); // expand the graph until dimension dim_max

	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", st.num_simplices());
	}

	// Sort the simplices in the order of the filtration
	st.initialize_filtration();

	// Compute the persistence diagram of the complex
	int p = 2; //characteristic of the coefficient field for homology
	double min_persistence = 0; //minimal length for persistent intervals
	computePersistenceGUDHI(st, p, min_persistence, maxdimension, persDgm,
			printProgress);

	// Output persistent diagram
	return Rcpp::List::create(
			concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
			concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
			StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// RipsDiag
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  dist           "euclidean" for Euclidean distance,
  *                            "arbitrary" for an arbitrary distance
  * @param[in]  library        Either "GUDHI", "Dionysus", or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
  *                            diagram. Diagram must point to enough memory space for
  *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
  *                            write nothing after.
  */
// [[Rcpp::export]]
Rcpp::List
RipsDiag(const Rcpp::NumericMatrix & X
       , const int                   maxdimension
       , const double                maxscale
       , const std::string         & dist
       , const std::string         & library
       , const bool                  location
       , const bool                  printProgress
	) {
	if (library[0] == 'G') {
		return RipsDiagGUDHI(X, maxdimension, maxscale, printProgress);
	}
	if (dist[0] == 'e') {
		return RipsDiagL2DionysusPhat(X, maxdimension, maxscale, library,
				location, printProgress);
	}
	else {
		return RipsDiagArbitDionysusPhat(X, maxdimension, maxscale, library,
				location, printProgress);
	}
}



//---------------------------------------------------------------------------------------------------------------------
// gudhi type definition
typedef Gudhi::Simplex_tree<>::Vertex_handle Simplex_tree_vertex;
typedef std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex > Alpha_shape_simplex_tree_map;
typedef std::pair<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex> Alpha_shape_simplex_tree_pair;
typedef std::vector< Simplex_tree_vertex > Simplex_tree_vector_vertex;

Vertex_list fromCell(const Cell_handle& ch)
{
	Vertex_list the_list;
	for (auto i = 0; i < 4; i++)
	{
		the_list.push_back(ch->vertex(i));
	}
	return the_list;
}
Vertex_list fromFacet(const Facet& fct)
{
	Vertex_list the_list;
	for (auto i = 0; i < 4; i++)
	{
		if (fct.second != i)
		{
			the_list.push_back(fct.first->vertex(i));
		}
	}
	return the_list;
}
Vertex_list fromEdge(const Edge_3& edg)
{
	Vertex_list the_list;
	for (auto i = 0; i < 4; i++)
	{
		if ((edg.second == i) || (edg.third == i))
		{
			the_list.push_back(edg.first->vertex(i));
		}
	}
	return the_list;
}
Vertex_list fromVertex(const Alpha_shape_3::Vertex_handle& vh)
{
	Vertex_list the_list;
	the_list.push_back(vh);
	return the_list;
}



// AlphaShapeDiag in GUDHI
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List
AlphaShapeDiagGUDHI(const Rcpp::NumericMatrix & X          //points to some memory space
                  , const bool                  printProgress
	) {
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	  int coeff_field_characteristic = 2;

	  float min_persistence = 0.0;



	  // Turn the input points into a range of points
	  std::list<Point_3> lp = RcppToCGALPoint3< std::list<Point_3> >(X);


	  // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
	  Alpha_shape_3 as(lp.begin(),lp.end(),0,Alpha_shape_3::GENERAL);
	  //std::cout << "Alpha shape computed in GENERAL mode" << std::endl;

	  // filtration with alpha values from alpha shape
	  std::vector<Object> the_objects;
	  std::vector<Alpha_value_type> the_alpha_values;

	  Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>( std::back_inserter(the_objects), std::back_inserter(the_alpha_values));

	  as.filtration_with_alpha_values(disp);
	  //std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;

	  Alpha_shape_3::size_type count_vertices = 0;
	  Alpha_shape_3::size_type count_edges    = 0;
	  Alpha_shape_3::size_type count_facets   = 0;
	  Alpha_shape_3::size_type count_cells    = 0;

	  // Loop on objects vector
	  Vertex_list vertex_list;
	  Gudhi::Simplex_tree<> simplex_tree;
	  Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
	  std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
	  int dim_max=0;
	  Filtration_value filtration_max=0.0;
	  for(auto object_iterator: the_objects)
	  {
	    // Retrieve Alpha shape vertex list from object
	    if (const Cell_handle* cell = CGAL::object_cast<Cell_handle>(&object_iterator))
	    {
	      vertex_list = fromCell(*cell);
	      count_cells++;
	      if (dim_max < 3) {
	        dim_max=3; // Cell is of dim 3
	      }
	    }
	    else if (const Facet* facet = CGAL::object_cast<Facet>(&object_iterator))
	    {
	      vertex_list = fromFacet(*facet);
	      count_facets++;
	      if (dim_max < 2) {
	        dim_max=2; // Facet is of dim 2
	      }
	    }
	    else if (const Edge_3* edge = CGAL::object_cast<Edge_3>(&object_iterator))
	    {
	      vertex_list = fromEdge(*edge);
	      count_edges++;
	      if (dim_max < 1) {
	        dim_max=1; // Edge_3 is of dim 1
	      }
	    }
	    else if (const Alpha_shape_3::Vertex_handle* vertex = CGAL::object_cast<Alpha_shape_3::Vertex_handle>(&object_iterator))
	    {
	      count_vertices++;
	      vertex_list = fromVertex(*vertex);
	    }
	    // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
	    Simplex_tree_vector_vertex the_simplex_tree;
	    for (auto the_alpha_shape_vertex:vertex_list)
	    {
	      Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
	      if (the_map_iterator == map_cgal_simplex_tree.end())
	      {
	        // alpha shape not found
	        Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();
	        //std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] not found - insert " << vertex << std::endl;
	        the_simplex_tree.push_back(vertex);
	        map_cgal_simplex_tree.insert(Alpha_shape_simplex_tree_pair(the_alpha_shape_vertex,vertex));
	      } else
	      {
	        // alpha shape found
	        Simplex_tree_vertex vertex = the_map_iterator->second;
	        //std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] found in " << vertex << std::endl;
	        the_simplex_tree.push_back(vertex);
	      }
	    }
	    // Construction of the simplex_tree
	    Filtration_value filtr = std::sqrt(*the_alpha_value_iterator);
	    //std::cout << "filtration = " << filtr << std::endl;
	    if (filtr > filtration_max) {
	      filtration_max = filtr;
	    }
	    simplex_tree.insert_simplex(the_simplex_tree, filtr);
		if (the_alpha_value_iterator != the_alpha_values.end()) {
		  ++the_alpha_value_iterator;
		}
		else {
		  //std::cout << "This shall not happen" << std::endl;
		}
	  }
	  simplex_tree.set_filtration(filtration_max);
	  simplex_tree.set_dimension(dim_max);

	  if (printProgress) {
		  Rprintf("# Generated complex of size: %d \n", simplex_tree.num_simplices());
	  }
	  //std::cout << "vertices \t\t" << count_vertices << std::endl;
	  //std::cout << "edges \t\t"    << count_edges << std::endl;
	  //std::cout << "facets \t\t"   << count_facets << std::endl;
	  //std::cout << "cells \t\t"    << count_cells << std::endl;


	  //std::cout << "Information of the Simplex Tree: " << std::endl;
	  //std::cout << "  Number of vertices = " << simplex_tree.num_vertices() << " ";
	  //std::cout << "  Number of simplices = " << simplex_tree.num_simplices() << std::endl << std::endl;
	  //std::cout << "  Dimension = " << simplex_tree.dimension() << " ";
	  //std::cout << "  filtration = " << simplex_tree.filtration() << std::endl << std::endl;
	  //std::cout << "Iterator on vertices: " << std::endl;
	  //for( auto vertex : simplex_tree.complex_vertex_range() )
	  //{ std::cout << vertex << " "; }

	  // Sort the simplices in the order of the filtration
	  simplex_tree.initialize_filtration();

	  //std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;

	// Compute the persistence diagram of the complex
	computePersistenceGUDHI(simplex_tree, coeff_field_characteristic,
			min_persistence, 2, persDgm, printProgress);

	// Output persistent diagram
	return Rcpp::List::create(
			concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
			concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
			StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}

// AlphaComplexDiag in GUDHI
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  maxalphasquare Threshold for the Alpha complex,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List
AlphaComplexDiagGUDHI(const Rcpp::NumericMatrix & X             //points to some memory space
                    , const bool                  printProgress
	) {
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	int coeff_field_characteristic = 2;

	float min_persistence = 0.0;

        using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>;
        using Point = Kernel::Point_d;

	// Turn the input points into a range of points
	std::list<Point> lp = RcppToCGALPointD< std::list< Point > >(X);

	//Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_points(lp, maxalphasquare);
  Gudhi::alphacomplex::Alpha_complex<Kernel>
      alpha_complex_from_points(lp, std::numeric_limits<double>::infinity());

	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", alpha_complex_from_points.num_simplices());
	}

	// Sort the simplices in the order of the filtration
	alpha_complex_from_points.initialize_filtration();

	// Compute the persistence diagram of the complex
	computePersistenceGUDHI(alpha_complex_from_points, coeff_field_characteristic,
			min_persistence, 2, persDgm, printProgress);

	// Output persistent diagram
	return Rcpp::List::create(
			concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
			concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
			StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}
