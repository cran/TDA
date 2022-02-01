#ifndef __FILTRATIONUTILS_H__
#define __FILTRATIONUTILS_H__

#include <vector>
#include <string>
#include <algorithm>
#include <functional>

#include <utilities/timer.h>



template< typename VertexVector, typename RealVector >
void filtrationSort(
    std::vector< VertexVector > & cmplx, RealVector & values) {

  std::vector< std::pair< double, unsigned > > vipairs(cmplx.size());

  {
    typename RealVector::iterator iValue = values.begin();
    unsigned idx = 0;
    typename std::vector< std::pair< double, unsigned > >::iterator iPair =
          vipairs.begin();
    for (; iPair != vipairs.end(); ++iPair, ++iValue, ++idx) {
      *iPair = std::make_pair(*iValue, idx);
    }
  }

  std::sort(vipairs.begin(), vipairs.end());

  {
    std::vector< VertexVector > cmplxTemp(cmplx.begin(), cmplx.end());

    typename std::vector< VertexVector >::iterator iCmplx = cmplx.begin();
    typename RealVector::iterator iValue = values.begin();
    typename std::vector< std::pair< double, unsigned > >::iterator iPair =
        vipairs.begin();
    for (; iPair != vipairs.end(); ++iPair, ++iCmplx, ++iValue) {
      *iCmplx = cmplxTemp[iPair->second];
      *iValue = iPair->first;
    }
  }
}



// funFiltration
/** \brief Interface for R code, construct the persistence diagram from the
*         filtration of the function values.
*
* @param[out] Rcpp::List  A list
* @param[in]  FUNvalues   The inlut function values
* @param[in]  cmplx       The input simplicial complex
*/template< typename RealVector, typename VertexVector >
inline void funFiltration(
    const RealVector            & FUNvalues,
    std::vector< VertexVector > & cmplx,
    std::vector< double >       & values
) {

  const unsigned nCmplx = cmplx.size();
  values = std::vector< double >(nCmplx);

  std::vector< double >::iterator iValue = values.begin();
  for (typename std::vector< VertexVector >::iterator iCmplx = cmplx.begin();
    iCmplx != cmplx.end(); ++iCmplx, ++iValue) {
    VertexVector cmplxVec(*iCmplx);
    typename VertexVector::iterator iCmplxVec = cmplxVec.begin();
    *iValue = FUNvalues[*iCmplxVec];
    for (; iCmplxVec != cmplxVec.end(); ++iCmplxVec) {
      *iValue = std::max(*iValue, FUNvalues[*iCmplxVec]);
    }
  }

  // sort
  filtrationSort(cmplx, values);
}



# endif // __FILTRATIONUTILS_H__

