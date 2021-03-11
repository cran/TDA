#include <string>
#include <vector>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>

#include <tdautils/gridUtils.h>



// AlphaShapeDiag
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *        shape complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  printProgress  Is progress printed?
 */
template< typename RealMatrix, typename Print >
void alphaShapeDiag(
  const RealMatrix  & X,             //points to some memory space
  const unsigned      nSample,
  const unsigned      nDim,
  const int           maxdimension,
  const std::string & libraryDiag,
  const bool          location,
  const bool          printProgress,
  const Print       & print,
  std::vector< std::vector< std::vector< double > > > & persDgm,
  std::vector< std::vector< std::vector< unsigned > > > & persLoc,
  std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle,
  RealMatrix        & coordinates
) {

  int coeff_field_characteristic = 2;

  float min_persistence = 0.0;

  // 2021-02-08, Jisu KIM
  // fixing [-Wclass-memaccess] warning
  //Gudhi::Simplex_tree<> smplxTree =
  //  AlphaShapeFiltrationGudhi< Gudhi::Simplex_tree<> >(
  //    X, printProgress, print, coordinates);
  Gudhi::Simplex_tree<> smplxTree(
    AlphaShapeFiltrationGudhi< Gudhi::Simplex_tree<> >(
      X, printProgress, print, coordinates));

  // Compute the persistence diagram of the complex
  if (libraryDiag[0] == 'G') {
    FiltrationDiagGudhi(
        smplxTree, coeff_field_characteristic, min_persistence, 2,
        printProgress, persDgm);
  }
  else if (libraryDiag[0] == 'D') {
    Fltr filtration = filtrationGudhiToDionysus< Fltr >(smplxTree);
    FiltrationDiagDionysus< Persistence >(
        filtration, maxdimension, location, printProgress, persDgm, persLoc,
        persCycle);
  }
  else {
    std::vector< phat::column > cmplx;
    std::vector< double > values;
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    filtrationGudhiToPhat< phat::column, phat::dimension >(
        smplxTree, cmplx, values, boundary_matrix);
    FiltrationDiagPhat(
        cmplx, values, boundary_matrix, maxdimension, location,
        printProgress, persDgm, persLoc, persCycle);
  }
}
