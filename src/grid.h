#include <cmath>
#include <vector>
#include <string>

// for grid
#include <tdautils/gridUtils.h>

#include <tdautils/filtrationUtils.h>
#include <tdautils/filtrationDiag.h>



template< typename RealVector, typename IntVector, typename RealMatrix >
inline void gridBy(
     const RealVector & lim, const double by, IntVector & dim,
     RealMatrix & grid, unsigned & nGrid) {

  unsigned nDim = (unsigned)(lim.size() / 2);

  dim = IntVector(nDim);
  for (unsigned iDim = 0; iDim < nDim; ++iDim) {
    dim[iDim] = (unsigned)std::floor((lim[2 * iDim + 1] - lim[2 * iDim]) / by) + 1;
  }

  nGrid = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies< double >());
  grid = RealMatrix(nGrid * nDim);

  RealVector gridVec(nDim);
  for (unsigned iDim = 0; iDim < nDim; ++iDim) {
    gridVec[iDim] = lim[2 * iDim];
  }

  for (unsigned iGrid = 0; iGrid < nGrid; ++iGrid) {
    for (unsigned iDim = 0; iDim < nDim; ++iDim) {
      grid[iGrid + iDim * nGrid] = gridVec[iDim];
    }

    gridVec[0] += by;
    for (unsigned iDim = 0; iDim < nDim; ++iDim) {
      if (gridVec[iDim] > lim[2 * iDim + 1] && iDim + 1 != nDim) {
        gridVec[iDim] = lim[2 * iDim];
        gridVec[iDim + 1] += by;
      }
    }
  }
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
template< typename RealVector, typename DimensionVector, typename Print,
          typename VertexVector >
inline void gridFiltration(
    const RealVector            & FUNvalues,
    const DimensionVector       & gridDim,
    const int                     maxdimension,
    const std::string           & decomposition,
    const bool                    printProgress,
    const Print                 & print,
    std::vector< VertexVector > & cmplx,
    std::vector< double >       & values
) {
#ifdef LOGGING
  //rlog::RLogInit(argc, argv);

  stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
  //stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
  //stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
#endif

  // Generate simplicial complex from function values and grid
  if (decomposition[0] == '5') {
    simplicesFromGrid(gridDim, maxdimension + 1, cmplx);
  }
  if (decomposition[0] == 'b') {
    simplicesFromGridBarycenter(gridDim, maxdimension + 1, cmplx);
  }
  if (printProgress) {
    print("# Generated complex of size: %d \n", cmplx.size());
  }

  // Sort the simplices with respect to function values
  //filtration.sort(Smplx::DataComparison());
  funFiltration(FUNvalues, cmplx, values);
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
template< typename VertexVector, typename RealVector, typename DimensionVector, typename Print >
inline void gridDiag(
    const RealVector      & FUNvalues,
    const DimensionVector & gridDim,
    const int               maxdimension,
    const std::string     & decomposition,
    const std::string     & library,
    const bool              location,
    const bool              printProgress,
    const Print           & print,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  //Fltr filtration;
  std::vector< VertexVector > cmplx;
  std::vector< double > values;

  gridFiltration(
      FUNvalues, gridDim, maxdimension, decomposition, printProgress, print,
      cmplx, values);

  // Compute the persistence diagram of the complex
  filtrationDiagSorted< VertexVector >(
      cmplx, values, maxdimension, library, location, printProgress, 0,
      persDgm, persLoc, persCycle);
}
