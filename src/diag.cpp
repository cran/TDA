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

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>

#include <tdautils/filtrationUtils.h>
#include <tdautils/filtrationDiag.h>

#include <alphaComplex.h>
#include <alphaShape.h>
#include <dtm.h>
#include <grid.h>
#include <kde.h>
#include <rips.h>



// GridFiltration by Brittany T. Fasy
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
Rcpp::List GridFiltration(
    const Rcpp::NumericVector & FUNvalues,
    const Rcpp::IntegerVector & gridDim,
    const int                   maxdimension,
    const std::string         & decomposition,
    const bool                  printProgress
) {
#ifdef LOGGING
  //rlog::RLogInit(argc, argv);

  stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
  //stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
  //stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
#endif

  std::vector< std::vector< unsigned > > cmplx;
  std::vector< double > values;

  gridFiltration(
      FUNvalues, gridDim, maxdimension, decomposition, printProgress, Rprintf,
      cmplx, values);

  return Rcpp::List::create(
      StlCmplxToRcpp< Rcpp::IntegerVector, Rcpp::List >(cmplx, 1),
      Rcpp::NumericVector(values.begin(), values.end()));
}



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
Rcpp::List GridDiag(
    const Rcpp::NumericVector & FUNvalues,
    const Rcpp::IntegerVector & gridDim,
    const int                   maxdimension,
    const std::string         & decomposition,
    const std::string         & library,
    const bool                  location,
    const bool                  printProgress
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

  gridDiag< std::vector< unsigned > >(
      FUNvalues, gridDim, maxdimension, decomposition, library, location,
      printProgress, Rprintf, persDgm, persLoc, persCycle);

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
Rcpp::NumericVector Kde(
    const Rcpp::NumericMatrix & X,
    const Rcpp::NumericMatrix & Grid,
    const double                h,
    const std::string         & kertype,
    const Rcpp::NumericVector & weight,
    const bool                  printProgress
	) {

  return kde(
      X, Grid, X.nrow(), Grid.ncol(), Grid.nrow(), h, kertype, weight,
      printProgress, Rprintf);
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

  return kdeDist(
      X, Grid, X.nrow(), Grid.ncol(), Grid.nrow(), h, weight, printProgress,
      Rprintf);
}



// distance to measure function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
Dtm(const Rcpp::NumericMatrix & knnDistance
  , const double                weightBound
  , const double                r
	) {

  return dtm< Rcpp::NumericVector >(
      knnDistance, knnDistance.nrow(), weightBound, r);
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

  return dtmWeight(
      knnDistance, knnDistance.nrow(), weightBound, r, knnIndex, weight);
}



// FiltrationDiag
/** \brief Interface for R code, construct the persistence diagram from the
 *         filtration.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  filtration     The input filtration
 * @param[in]  maxdimension   Max dimension of the homological features to be
 *                            computed.
 * @param[in]  library        Either "GUDHI", "Dionysus", or "PHAT"
 * @param[in]  location       Are location of birth point, death point, and
 *                            representative cycles returned?
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List FiltrationDiag(
    const Rcpp::List  & filtration,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress
) {

  Rcpp::List cmplx(filtration[0]);
  Rcpp::NumericVector values(filtration[1]);
  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;

  filtrationDiag< Rcpp::IntegerVector >(
      cmplx, values, maxdimension, library, location, printProgress, 1,
      persDgm, persLoc, persCycle);

  // Output persistent diagram
  return Rcpp::List::create(
    concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
    concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
    StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// funFiltration
/** \brief Interface for R code, construct the persistence diagram from the
 *         filtration of the function values.
 *
 * @param[out] Rcpp::List  A list
 * @param[in]  FUNvalues   The inlut function values
 * @param[in]  cmplx       The input simplicial complex
 */
// [[Rcpp::export]]
Rcpp::List FunFiltration(
    const Rcpp::NumericVector & FUNvalues,
    const Rcpp::List          & cmplx
) {

  std::vector< std::vector< unsigned > > funCmplx =
      RcppCmplxToStl< std::vector< unsigned >, Rcpp::IntegerVector >(cmplx, 1);
  std::vector< double > values;

  funFiltration(FUNvalues, funCmplx, values);

  return Rcpp::List::create(
      StlCmplxToRcpp< Rcpp::IntegerVector, Rcpp::List >(funCmplx, 1),
      Rcpp::NumericVector(values.begin(), values.end()));
}



// RipsFiltration
/** \brief Interface for R code, construct the rips filtration on the input
 *         set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              Either an nxd matrix of coordinates,
 *                            or an nxn matrix of distances of points
 * @param[in]  maxdimension   Max dimension of the homological features to be computed.
 * @param[in]  maxscale       Threshold for the Rips complex
 * @param[in]  dist           "euclidean" for Euclidean distance,
 *                            "arbitrary" for an arbitrary distance
 * @param[in]  library        Either "GUDHI" or "Dionysus"
 * @param[in]  printProgress  Is progress printed?
 * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
 *                            diagram. Diagram must point to enough memory space for
 *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
 *                            write nothing after.
 */
// [[Rcpp::export]]
Rcpp::List RipsFiltration(
    const Rcpp::NumericMatrix & X,
    const int                   maxdimension,
    const double                maxscale,
    const std::string         & dist,
    const std::string         & library,
    const bool                  printProgress
) {

  Rcpp::List cmplx;
  Rcpp::NumericVector values;
  Rcpp::List boundary;

  ripsFiltration< Rcpp::IntegerVector >(
      X, X.nrow(), X.ncol(), maxdimension, maxscale, dist, library,
      printProgress, Rprintf, cmplx, values, boundary);

  return Rcpp::List::create(cmplx, values, boundary);
}



// RipsDiag
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  dist           "euclidean" for Euclidean distance,
  *                            "arbitrary" for an arbitrary distance
  * @param[in]  libraryFiltration  Either "GUDHI" or "Dionysus"
  * @param[in]  libraryDiag        Either "GUDHI", "Dionysus", or "PHAT"
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
       , const std::string         & libraryFiltration
       , const std::string         & libraryDiag
       , const bool                  location
       , const bool                  printProgress
) {

	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

  ripsDiag(
      X, X.nrow(), X.ncol(), maxdimension, maxscale, dist, libraryFiltration,
      libraryDiag, location, printProgress, Rprintf, persDgm, persLoc,
      persCycle);

	// Output persistent diagram
	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
		concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
		StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// AlphaShapeFiltration in GUDHI
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *         shape complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List AlphaShapeFiltration(
  const Rcpp::NumericMatrix & X,          //points to some memory space
  const bool                  printProgress
) {

  Rcpp::NumericMatrix coordinates;

  Gudhi::Simplex_tree<> smplxTree =
      AlphaShapeFiltrationGudhi< Gudhi::Simplex_tree<> >(
          X, printProgress, Rprintf, coordinates);
  Rcpp::List filtration =
      filtrationGudhiToRcpp< Rcpp::List, Rcpp::NumericVector >(smplxTree);
  filtration.push_back(coordinates);
  return filtration;
}



// AlphaShapeDiag
/** \brief Interface for R code, construct the persistence diagram of the alpha
  *        shape complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List AlphaShapeDiag(
    const Rcpp::NumericMatrix & X,             //points to some memory space
    const int                   maxdimension,
    const std::string         & libraryDiag,
    const bool                  location,
    const bool                  printProgress
	) {

  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;
  Rcpp::NumericMatrix coordinates;

  alphaShapeDiag(
      X, X.nrow(), X.ncol(), maxdimension, libraryDiag, location,
      printProgress, Rprintf, persDgm, persLoc, persCycle, coordinates);

  // Output persistent diagram
  return Rcpp::List::create(
      concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
      concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
      StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle),
      coordinates);
}



// AlphaComplexFiltration
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *         complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  maxalphasquare Threshold for the Alpha complex,
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List AlphaComplexFiltration(
  const Rcpp::NumericMatrix & X,             //points to some memory space
  const bool                  printProgress
) {

  Gudhi::Simplex_tree<> smplxTree =
      AlphaComplexFiltrationGudhi< Gudhi::Simplex_tree<> >(
          X, printProgress, Rprintf);
  // 2018-08-04
  // switching back to original code
  return filtrationGudhiToRcpp< Rcpp::List, Rcpp::NumericVector >(smplxTree);
}



// AlphaComplexDiag
/** \brief Interface for R code, construct the persistence diagram of the alpha
  *        complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  maxalphasquare Threshold for the Alpha complex,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List AlphaComplexDiag(
    const Rcpp::NumericMatrix & X,             //points to some memory space
    const int                   maxdimension,
    const std::string         & libraryDiag,
    const bool                  location,
    const bool                  printProgress
	) {

  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;

  alphaComplexDiag(
    X, X.nrow(), X.ncol(), maxdimension, libraryDiag, location, printProgress,
    Rprintf, persDgm, persLoc, persCycle);

  // Output persistent diagram
  return Rcpp::List::create(
      concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
      concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
      StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}
