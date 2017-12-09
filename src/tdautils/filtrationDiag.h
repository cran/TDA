#ifndef __FILTRATIONDIAG_H__
#define __FILTRATIONDIAG_H__

#include <vector>
#include <string>
#include <algorithm>

#include <tdautils/filtrationUtils.h>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>

// for grid
#include <tdautils/gridUtils.h>

#include <iostream>



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
// TODO: see whether IntegerVector in template is deducible
template< typename VertexVector, typename VectorList, typename RealVector >
inline void filtrationDiagSorted(
    VectorList        & cmplx,
    RealVector        & values,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress,
    const unsigned      idxShift,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  if (library[0] == 'G') {
    int coeff_field_characteristic = 2;
    double min_persistence = 0.0;
    Gudhi::Simplex_tree<> smplxTree = filtrationTdaToGudhi<
        VertexVector, Gudhi::Simplex_tree<> >(
            cmplx, values, idxShift);
    FiltrationDiagGudhi(
        smplxTree, coeff_field_characteristic, min_persistence, maxdimension,
        printProgress, persDgm);
  }
  else if (library[0] == 'D') {
    FiltrationDiagDionysus< Persistence >(
        filtrationTdaToDionysus< VertexVector, Fltr >(
            cmplx, values, idxShift),
        maxdimension, location, printProgress, persDgm, persLoc, persCycle);
  }
  else {

    std::vector< phat::column > cmplxPhat(cmplx.size());
    typename VectorList::iterator iCmplx = cmplx.begin();
    std::vector< phat::column >::iterator iPhat = cmplxPhat.begin();
    for (; iCmplx != cmplx.end(); ++iCmplx, ++iPhat) {
      VertexVector cmplxVec(*iCmplx);
      *iPhat = phat::column(cmplxVec.begin(), cmplxVec.end());
    }

    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    filtrationDionysusToPhat< phat::column, phat::dimension >(
        filtrationTdaToDionysus< phat::column, Fltr >(
            cmplxPhat, values, idxShift),
        cmplxPhat, values, boundary_matrix);
    FiltrationDiagPhat(cmplxPhat, values, boundary_matrix,
        maxdimension, location, printProgress, persDgm, persLoc, persCycle);
  }
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
// TODO: see whether IntegerVector in template is deducible
template< typename VertexVector, typename VectorList, typename RealVector >
inline void filtrationDiag(
    VectorList        & cmplx,
    RealVector        & values,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress,
    const unsigned      idxShift,
    std::vector< std::vector< std::vector< double > > > & persDgm,
    std::vector< std::vector< std::vector< unsigned > > > & persLoc,
    std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {

  if (std::is_sorted(values.begin(), values.end())) {
    filtrationDiagSorted< VertexVector >(
        cmplx, values, maxdimension, library, location, printProgress,
        idxShift, persDgm, persLoc, persCycle);
  }
  else {
    std::vector< std::vector< unsigned > > cmplxTemp = 
        RcppCmplxToStl< std::vector< unsigned >, VertexVector >(cmplx, 0);
    std::vector< double > valuesTemp(values.begin(), values.end());
    filtrationSort(cmplxTemp, valuesTemp);
    filtrationDiagSorted< std::vector< unsigned > >(
        cmplxTemp, valuesTemp, maxdimension, library, location, printProgress,
        idxShift, persDgm, persLoc, persCycle);
  }
}



# endif // __FILTRATIONDIAG_H__
