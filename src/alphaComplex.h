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



// AlphaComplexDiag
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *        complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  maxalphasquare Threshold for the Alpha complex,
 * @param[in]  printProgress  Is progress printed?
 */
template< typename RealMatrix, typename Print >
void alphaComplexDiag(
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
  std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>;
  using Point = Kernel::Point_d;

  int coeff_field_characteristic = 2;

  float min_persistence = 0.0;

  Gudhi::Simplex_tree<> alphaCmplx =
      AlphaComplexFiltrationGudhi< Gudhi::Simplex_tree<> >(
          X, printProgress, print);

  // 2018-08-04
  // switching back to original code

  // Compute the persistence diagram of the complex
  if (libraryDiag[0] == 'G') {
    // 2018-08-04
    // switching back to original code
    FiltrationDiagGudhi(
        alphaCmplx, coeff_field_characteristic, min_persistence, 2,
        printProgress, persDgm);
  }
  else if (libraryDiag[0] == 'D') {
    // 2018-08-04
    // switching back to original code
    Fltr filtration = filtrationGudhiToDionysus< Fltr >(alphaCmplx);
    FiltrationDiagDionysus< Persistence >(
        filtration, maxdimension, location, printProgress, persDgm, persLoc,
        persCycle);
  }
  else {
    // 2018-08-04
    // switching back to original code
    std::vector< phat::column > cmplx;
    std::vector< double > values;
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    filtrationGudhiToPhat< phat::column, phat::dimension >(
        alphaCmplx, cmplx, values, boundary_matrix);
    FiltrationDiagPhat(
        cmplx, values, boundary_matrix, maxdimension, location,
        printProgress, persDgm, persLoc, persCycle);
  }
}
