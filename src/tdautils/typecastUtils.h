#ifndef __TYPECASTUTILS_H__
#define __TYPECASTUTILS_H__

#include <vector>
#include <map>
#include <algorithm>



template<typename PersistenceDiagram, typename RcppMatrix>
inline PersistenceDiagram RcppToDionysus(const RcppMatrix& rcppMatrix) {
	PersistenceDiagram dionysusDiagram;
	const unsigned rowNum = rcppMatrix.nrow();
	for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx)
	{
		dionysusDiagram.push_back(typename PersistenceDiagram::Point(
				rcppMatrix[rowIdx + 0 * rowNum], rcppMatrix[rowIdx + 1 * rowNum]));
	}
	return dionysusDiagram;
}



template< typename StlMatrix, typename RealMatrix >
inline StlMatrix TdaToStl(const RealMatrix & rcppMatrix,
    const unsigned nRow, const unsigned nCol, bool is_row_names = false) {

  if (is_row_names) {
    StlMatrix stlMatrix(nRow, typename StlMatrix::value_type(nCol + 1));
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      stlMatrix[rowIdx][0] = rowIdx + 1;
    }
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      for (unsigned colIdx = 0; colIdx < nCol; ++colIdx) {
        stlMatrix[rowIdx][colIdx + 1] = rcppMatrix[rowIdx + colIdx * nRow];
      }
    }
    return stlMatrix;
  }
  else {
    StlMatrix stlMatrix(nRow, typename StlMatrix::value_type(nCol));
    for (unsigned rowIdx = 0; rowIdx < nRow; ++rowIdx) {
      for (unsigned colIdx = 0; colIdx < nCol; ++colIdx) {
        stlMatrix[rowIdx][colIdx] = rcppMatrix[rowIdx + colIdx * nRow];
      }
    }
    return stlMatrix;
  }
}



template< typename StlMatrix, typename RcppMatrix >
inline StlMatrix RcppToStl(const RcppMatrix& rcppMatrix,
		bool is_row_names = false) {

	const unsigned rowNum = rcppMatrix.nrow();
	const unsigned colNum = rcppMatrix.ncol();
	if (is_row_names) {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum + 1));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			stlMatrix[rowIdx][0] = rowIdx + 1;
		}
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx + 1] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	else {
		StlMatrix stlMatrix(rowNum, typename StlMatrix::value_type(colNum));
		for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
			for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
				stlMatrix[rowIdx][colIdx] = rcppMatrix[rowIdx + colIdx * rowNum];
			}
		}
		return stlMatrix;
	}
	
}



template< typename CGALPoint3List, typename RcppMatrix >
inline CGALPoint3List RcppToCGALPoint3(const RcppMatrix& rcppMatrix) {

	const unsigned rowNum = rcppMatrix.nrow();
  CGALPoint3List cGALPoint3List;
	for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
    cGALPoint3List.push_back(
				typename CGALPoint3List::value_type(rcppMatrix[rowIdx],
					rcppMatrix[rowIdx + rowNum], rcppMatrix[rowIdx + 2 * rowNum]));
	}
	return cGALPoint3List;
}



template< typename CGALPointDList, typename RcppMatrix >
inline CGALPointDList RcppToCGALPointD(const RcppMatrix& rcppMatrix) {

  const unsigned rowNum = rcppMatrix.nrow();
  const unsigned colNum = rcppMatrix.ncol();
  CGALPointDList cGALPointDList;
  std::vector< double > pointD(colNum);

  for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
    for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
      pointD[colIdx] = rcppMatrix[rowIdx + colIdx * rowNum];
    }
    cGALPointDList.push_back(
        typename CGALPointDList::value_type(pointD.size(), pointD.begin(),
          pointD.end()));
  }
  return cGALPointDList;
}



template< typename VertexVector, typename RcppVector, typename RcppList >
inline std::vector< VertexVector > RcppCmplxToStl(
    const RcppList & rcppCmplx, const int idxShift) {

  const unsigned nCmplx = rcppCmplx.size();
  std::vector< VertexVector > stlCmplx(nCmplx);

  typename RcppList::const_iterator iRcppVec = rcppCmplx.begin();
  typename std::vector< VertexVector >::iterator iStlVec = stlCmplx.begin();
  for (; iRcppVec != rcppCmplx.end(); ++iRcppVec, ++iStlVec) {
    RcppVector cmplxVec(*iRcppVec);
    *iStlVec = VertexVector(cmplxVec.size());

    typename RcppVector::const_iterator iRcpp = cmplxVec.begin();
    typename VertexVector::iterator iStl = iStlVec->begin();
    for (; iRcpp != cmplxVec.end(); ++iRcpp, ++iStl) {
      *iStl = *iRcpp - idxShift;
    }
  }

  return stlCmplx;
}



template< typename RcppVector, typename RcppList, typename VectorList >
inline RcppList StlCmplxToRcpp(
    const VectorList & stlCmplx, const int idxShift) {

  const unsigned nCmplx = stlCmplx.size();
  RcppList rcppCmplx(nCmplx);

  typename VectorList::const_iterator iStlVec = stlCmplx.begin();
  typename RcppList::iterator iRcppVec = rcppCmplx.begin();
  for (; iStlVec != stlCmplx.end(); ++iStlVec, ++iRcppVec) {
    RcppVector cmplxVec(iStlVec->size());

    typename VectorList::value_type::const_iterator iStl = iStlVec->begin();
    typename RcppVector::iterator iRcpp = cmplxVec.begin();
    for (; iStl != iStlVec->end(); ++iStl, ++iRcpp) {
      *iRcpp = *iStl + idxShift;
    }
    *iRcppVec = cmplxVec;
  }

  return rcppCmplx;
}



template<typename RcppMatrix, typename StlMatrix>
inline RcppMatrix concatStlToRcpp(const std::vector< StlMatrix >& stlMatrices,
		bool includeIndex, unsigned colNum) {
	unsigned rowNum = 0;

	typename std::vector< StlMatrix >::const_iterator vecItr;
	for (vecItr = stlMatrices.begin(); vecItr != stlMatrices.end(); ++vecItr) {
		rowNum += vecItr->size();
	}
	RcppMatrix rcppMatrix(rowNum, colNum);

	unsigned vecIdx, rowIdx, colIdx;
	for (vecIdx = 0, rowIdx = 0; vecIdx < stlMatrices.size(); ++vecIdx) {
		typename StlMatrix::const_iterator matItr;
		for (matItr = stlMatrices[vecIdx].begin();
				matItr != stlMatrices[vecIdx].end(); ++matItr, ++rowIdx) {
			if (includeIndex) {
				rcppMatrix[rowIdx] = vecIdx;
				for (colIdx = 0; colIdx < colNum - 1; ++colIdx) {
					rcppMatrix[rowIdx + (colIdx + 1) * rowNum] = (*matItr)[colIdx];
				}
			}
			else {
				for (colIdx = 0; colIdx < colNum; ++colIdx) {
					rcppMatrix[rowIdx + colIdx * rowNum] = (*matItr)[colIdx];
				}
			}
		}
	}

	return rcppMatrix;
}



template<typename RcppList, typename RcppVector, typename StlSet>
inline RcppList StlToRcppList(
		const std::vector< std::vector< StlSet > >& stlSets) {
	unsigned rowNum = 0;

	typename std::vector< std::vector< StlSet > >::const_iterator vecsItr;
	for (vecsItr = stlSets.begin(); vecsItr != stlSets.end(); ++vecsItr) {
		rowNum += vecsItr->size();
	}
	RcppList rcppList(rowNum);

	typename RcppList::iterator listItr;
	typename std::vector< StlSet >::const_iterator setVecItr;
	typename StlSet::const_iterator setItr;
	unsigned setIdx;
	for (vecsItr = stlSets.begin(), listItr = rcppList.begin();
			vecsItr != stlSets.end(); ++vecsItr) {

		for (setVecItr = vecsItr->begin(); setVecItr != vecsItr->end(); ++setVecItr, ++listItr) {
			RcppVector rcppVec(setVecItr->size());
			for (setIdx = 0, setItr = setVecItr->begin(); setItr != setVecItr->end(); ++setItr, ++setIdx) {
				rcppVec[setIdx] = *setItr;
			}
			*listItr = rcppVec;
		}
	}

	return rcppList;
}



template< typename RcppList, typename RcppMatrix, typename StlVector >
inline RcppList StlToRcppMatrixList(
	const std::vector< std::vector< std::vector< StlVector > > >& stlArrays) {
	unsigned listNum = 0;

	typename std::vector< std::vector< std::vector< StlVector > > >::const_iterator vecsItr;
	for (vecsItr = stlArrays.begin(); vecsItr != stlArrays.end(); ++vecsItr) {
		listNum += vecsItr->size();
	}
	RcppList rcppList(listNum);

	typename RcppList::iterator listItr;
	typename std::vector< std::vector< StlVector > >::const_iterator matrixItr;
	typename std::vector< StlVector >::const_iterator rowItr;
	typename StlVector::const_iterator colItr;
	unsigned rowIdx, colIdx, rowNum;
	for (vecsItr = stlArrays.begin(), listItr = rcppList.begin();
	vecsItr != stlArrays.end(); ++vecsItr) {

		for (matrixItr = vecsItr->begin(); matrixItr != vecsItr->end();
		++matrixItr, ++listItr) {
			rowNum = matrixItr->size();
			if (rowNum != 0) {
				RcppMatrix rcppMatrix(rowNum, (*matrixItr)[0].size());
				for (rowIdx = 0, rowItr = matrixItr->begin();
				rowItr != matrixItr->end(); ++rowIdx, ++rowItr) {
					for (colIdx = 0, colItr = rowItr->begin(); colItr != rowItr->end();
					++colIdx, ++colItr) {
						rcppMatrix[rowIdx + colIdx * rowNum] = *colItr;
					}
				}
				*listItr = rcppMatrix;
			}
			else {
				RcppMatrix rcppMatrix(0, 0);
				*listItr = rcppMatrix;
			}
		}
	}

	return rcppList;
}



template< typename SimplexHandle, typename SimplexTree, typename RealVector >
void filtrationGudhiOne(
    // 2021-02-08, Jisu KIM
    // fixing [-Wclass-memaccess] warning
    //const SimplexHandle & sh, SimplexTree & smplxTree, const int idxShift,
    const SimplexHandle & sh, const SimplexTree & smplxTree, const int idxShift,
    // 2018-08-04
    // switching back to original code
    RealVector & cmplxVec, double & value, RealVector & boundaryVec) {

  const unsigned nVtx = smplxTree.dimension(sh) + 1;

  cmplxVec = RealVector(nVtx);
  const typename SimplexTree::Simplex_vertex_range & vtxRange =
      smplxTree.simplex_vertex_range(sh);
  typename RealVector::iterator iCmplxVec = cmplxVec.begin();
  for (typename SimplexTree::Simplex_vertex_iterator iVtx = vtxRange.begin();
       iVtx != vtxRange.end(); ++iVtx, ++iCmplxVec) {
    // R is 1-base, while C++ is 0-base
    *iCmplxVec = *iVtx + idxShift;
  }

  value = SimplexTree::filtration(sh);

  // 2018-08-04
  // switching back to original code

  // might need to change for cubical complex
  if (nVtx > 1) {
    boundaryVec = RealVector(nVtx);
  }
  const typename SimplexTree::Boundary_simplex_range & smplxRange =
      smplxTree.boundary_simplex_range(sh);

  typename RealVector::iterator iBdyVec = boundaryVec.begin();
  for (typename SimplexTree::Boundary_simplex_iterator iBdySpx =
      smplxRange.begin(); iBdySpx != smplxRange.end(); ++iBdySpx, ++iBdyVec) {
    // R is 1-base, while C++ is 0-base
    *iBdyVec = SimplexTree::key(*iBdySpx) + idxShift;
  }
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename IntegerVector, typename SimplexTree, typename VectorList,
          typename RealVector >
inline void filtrationGudhiToTda(
    // 2021-02-08, Jisu KIM
    // fixing [-Wclass-memaccess] warning
    //SimplexTree & smplxTree, VectorList & cmplx, RealVector & values,
    const SimplexTree & smplxTree, VectorList & cmplx, RealVector & values,
    // 2018-08-04
    // switching back to original code
    VectorList & boundary) {

  const unsigned nFltr = smplxTree.num_simplices();
  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary = VectorList(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();
  typename VectorList::iterator iBdy = boundary.begin();

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
    smplxTree.filtration_simplex_range();

  unsigned iFill = 0;
  for (typename SimplexTree::Filtration_simplex_iterator iSt =
    fltrGudhi.begin(); iSt != fltrGudhi.end();
    ++iSt, ++iCmplx, ++iValue, ++iBdy) {

    // Below two lines are only needed for computing boundary
    // 2021-02-08, Jisu KIM
    // temporarily fixing [-Wclass-memaccess] warning
    //smplxTree.assign_key(*iSt, iFill);
    //iFill++;

    IntegerVector cmplxVec;
    IntegerVector boundaryVec;
    // 2018-08-04
    // switching back to original code
    filtrationGudhiOne(*iSt, smplxTree, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;
  }
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename RcppList, typename RcppVector, typename SimplexTree >
// 2018-08-04
// switching back to original code
// 2021-02-08, Jisu KIM
// fixing [-Wclass-memaccess] warning
//inline RcppList filtrationGudhiToRcpp(SimplexTree & smplxTree) {
inline RcppList filtrationGudhiToRcpp(const SimplexTree & smplxTree) {

  const unsigned nFltr = smplxTree.num_simplices();

  RcppList cmplx(nFltr);
  RcppVector values(nFltr);
  RcppList boundary(nFltr);
  typename RcppList::iterator iCmplx = cmplx.begin();
  typename RcppVector::iterator iValue = values.begin();
  typename RcppList::iterator iBdy = boundary.begin();

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
      smplxTree.filtration_simplex_range();

  unsigned iFill = 0;
  for (typename SimplexTree::Filtration_simplex_iterator iSt =
       fltrGudhi.begin(); iSt != fltrGudhi.end();
       ++iSt, ++iCmplx, ++iValue, ++iBdy) {

    // Below two lines are only needed for computing boundary
    // 2021-02-08, Jisu KIM
    // temporarily fixing [-Wclass-memaccess] warning
    //smplxTree.assign_key(*iSt, iFill);
    //iFill++;

    RcppVector cmplxVec;
    RcppVector boundaryVec;
    // 2018-08-04
    // switching back to original code
    filtrationGudhiOne(*iSt, smplxTree, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;
  }

  return RcppList::create(cmplx, values, boundary);
}



template< typename SimplexTree, typename RcppVector, typename RcppList >
inline SimplexTree filtrationRcppToGudhi(const RcppList & rcppList) {

  const RcppList rcppComplex(rcppList[0]);
  const RcppVector rcppValue(rcppList[1]);
  SimplexTree smplxTree;

  typename RcppList::const_iterator iCmplx = rcppComplex.begin();
  typename RcppVector::const_iterator iValue = rcppValue.begin();
  for (; iCmplx != rcppComplex.end(); ++iCmplx, ++iValue) {
    const RcppVector rcppVec(*iCmplx);
    RcppVector gudhiVec(rcppVec.size());
    typename RcppVector::const_iterator iRcpp = rcppVec.begin();
    typename RcppVector::iterator iGudhi = gudhiVec.begin();
    for (; iRcpp != rcppVec.end(); ++iRcpp, ++iGudhi) {
      // R is 1-base, while C++ is 0-base
      *iGudhi = *iRcpp - 1;
    }
    smplxTree.insert_simplex(gudhiVec, *iValue);
  }

  return smplxTree;
}



template< typename IntegerVector, typename SimplexTree, typename VectorList,
          typename RealVector >
inline SimplexTree filtrationTdaToGudhi(
    const VectorList & cmplx, const RealVector & values, 
    const unsigned idxShift) {

  SimplexTree smplxTree;

  typename VectorList::const_iterator iCmplx = cmplx.begin();
  typename RealVector::const_iterator iValue = values.begin();
  for (; iCmplx != cmplx.end(); ++iCmplx, ++iValue) {
    const IntegerVector tdaVec(*iCmplx);
    IntegerVector gudhiVec(tdaVec.size());
    typename IntegerVector::const_iterator iTda = tdaVec.begin();
    typename IntegerVector::iterator iGudhi = gudhiVec.begin();
    for (; iTda != tdaVec.end(); ++iTda, ++iGudhi) {
      // R is 1-base, while C++ is 0-base
      *iGudhi = *iTda - idxShift;
    }
    smplxTree.insert_simplex(gudhiVec, *iValue);
  }

  return smplxTree;
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename Filtration, typename SimplexTree >
inline Filtration filtrationGudhiToDionysus(SimplexTree & smplxTree) {

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
      smplxTree.filtration_simplex_range();
  Filtration fltrDionysus;
  unsigned iFill = 0;

  for (typename SimplexTree::Filtration_simplex_iterator iSt =
       fltrGudhi.begin(); iSt != fltrGudhi.end(); ++iSt) {

    // Below two lines are only needed for computing boundary
    smplxTree.assign_key(*iSt, iFill);
    iFill++;

    std::vector< double > cmplxVec;
    double value;
    std::vector< double > boundaryVec;
    filtrationGudhiOne(*iSt, smplxTree, 0, cmplxVec, value, boundaryVec);

    fltrDionysus.push_back(typename Filtration::Simplex(
      cmplxVec.begin(), cmplxVec.end(), value));
  }

  return fltrDionysus;
}



// TODO : see whether 'const SimplexTree &' is possible
template< typename Column, typename Dimension, typename SimplexTree,
          typename VectorList, typename RealVector, typename Boundary >
inline void filtrationGudhiToPhat(
    SimplexTree & smplxTree, VectorList & cmplx, RealVector & values,
    Boundary & boundary_matrix) {

  const unsigned nFltr = smplxTree.num_simplices();

  const typename SimplexTree::Filtration_simplex_range & fltrGudhi =
      smplxTree.filtration_simplex_range();
  unsigned iFill = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary_matrix.set_num_cols(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();
  
  unsigned iCol = 0;
  for (typename SimplexTree::Filtration_simplex_iterator iSt =
       fltrGudhi.begin(); iSt != fltrGudhi.end();
       ++iSt, ++iCmplx, ++iValue, ++iCol) {

    // Below two lines are only needed for computing boundary
    smplxTree.assign_key(*iSt, iFill);
    iFill++;

    Column cmplxVec;
    Column boundary_indices;
    filtrationGudhiOne(
        *iSt, smplxTree, 0, cmplxVec, *iValue, boundary_indices);
    *iCmplx = cmplxVec;

    std::sort(boundary_indices.begin(), boundary_indices.end());
    boundary_matrix.set_col(iCol, boundary_indices);
    Dimension dim_of_column = smplxTree.dimension(*iSt);
    boundary_matrix.set_dim(iCol, dim_of_column);
  }
}



template< typename Simplex, typename SimplexMap, typename RealVector >
inline void filtrationDionysusOne(
  const Simplex & c, const SimplexMap & simplex_map, const int idxShift,
  RealVector & cmplxVec, double & value, RealVector & boundaryVec) {

  const unsigned nVtx = c.dimension() + 1;

  cmplxVec = RealVector(nVtx);
  typename RealVector::iterator iCmplxVec = cmplxVec.begin();
  for (typename Simplex::VertexContainer::const_iterator vit =
       c.vertices().begin(); vit != c.vertices().end(); ++vit, ++iCmplxVec) {
    // R is 1-base, while C++ is 0-base
    *iCmplxVec = *vit + idxShift;
  }

  value = c.data();

  // might need to change for cubical complex
  if (nVtx > 1) {
    boundaryVec = RealVector(nVtx);
  }
  typename RealVector::iterator iBdyVec = boundaryVec.begin();
  for (typename Simplex::BoundaryIterator bit = c.boundary_begin();
       bit != c.boundary_end(); ++bit, ++iBdyVec) {
    // R is 1-base, while C++ is 0-base
    *iBdyVec = simplex_map.find(*bit)->second + idxShift;
  }
}



template< typename IntegerVector, typename Filtration, typename VectorList,
          typename RealVector >
inline void filtrationDionysusToTda(
    const Filtration & filtration, VectorList & cmplx, RealVector & values,
    VectorList & boundary) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, unsigned,
      typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary = VectorList(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();
  typename VectorList::iterator iBdy = boundary.begin();

  for (typename Filtration::Index it = filtration.begin();
      it != filtration.end(); ++it, ++iCmplx, ++iValue, ++iBdy) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    IntegerVector cmplxVec;
    IntegerVector boundaryVec;
    filtrationDionysusOne(c, simplex_map, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;

    simplex_map.insert(typename
        std::map< typename Filtration::Simplex, unsigned >::value_type(
        c, size_of_simplex_map++));
  }
}



template< typename RcppList, typename RcppVector, typename Filtration >
inline RcppList filtrationDionysusToRcpp(const Filtration & filtration) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, unsigned,
    typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;

  RcppList cmplx(nFltr);
  RcppVector values(nFltr);
  RcppList boundary(nFltr);
  typename RcppList::iterator iCmplx = cmplx.begin();
  typename RcppVector::iterator iValue = values.begin();
  typename RcppList::iterator iBdy = boundary.begin();

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it, ++iCmplx, ++iValue, ++iBdy) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    RcppVector cmplxVec;
    RcppVector boundaryVec;
    filtrationDionysusOne(c, simplex_map, 1, cmplxVec, *iValue, boundaryVec);
    *iCmplx = cmplxVec;
    *iBdy = boundaryVec;

    simplex_map.insert(typename
        std::map< typename Filtration::Simplex, unsigned >::value_type(
            c, size_of_simplex_map++));
  }

  return RcppList::create(cmplx, values, boundary);
}



template< typename IntegerVector, typename Filtration, typename VectorList,
          typename RealVector >
inline Filtration filtrationTdaToDionysus(
    const VectorList & cmplx, const RealVector & values,
    const unsigned idxShift) {

  Filtration filtration;

  typename VectorList::const_iterator iCmplx = cmplx.begin();
  typename RealVector::const_iterator iValue = values.begin();
  for (; iCmplx != cmplx.end(); ++iCmplx, ++iValue) {
    const IntegerVector tdaVec(*iCmplx);
    IntegerVector dionysusVec(tdaVec.size());
    typename IntegerVector::const_iterator iTda = tdaVec.begin();
    typename IntegerVector::iterator iDionysus = dionysusVec.begin();
    for (; iTda != tdaVec.end(); ++iTda, ++iDionysus) {
      // R is 1-base, while C++ is 0-base
      *iDionysus = *iTda - idxShift;
    }
    filtration.push_back(typename Filtration::Simplex(
        dionysusVec.begin(), dionysusVec.end(), *iValue));
  }

  return filtration;
}



template< typename Filtration, typename RcppVector, typename RcppList >
inline Filtration filtrationRcppToDionysus(const RcppList & rcppList) {

  const RcppList rcppComplex(rcppList[0]);
  const RcppVector rcppValue(rcppList[1]);
  Filtration filtration;

  typename RcppList::const_iterator iCmplx = rcppComplex.begin();
  typename RcppVector::const_iterator iValue = rcppValue.begin();
  for (; iCmplx != rcppComplex.end(); ++iCmplx, ++iValue) {
    const RcppVector rcppVec(*iCmplx);
    RcppVector dionysusVec(rcppVec.size());
    typename RcppVector::const_iterator iRcpp = rcppVec.begin();
    typename RcppVector::iterator iDionysus = dionysusVec.begin();
    for (; iRcpp != rcppVec.end(); ++iRcpp, ++iDionysus) {
      // R is 1-base, while C++ is 0-base
      *iDionysus = *iRcpp - 1;
    }
    filtration.push_back(typename Filtration::Simplex(
        dionysusVec.begin(), dionysusVec.end(), *iValue));
  }

  return filtration;
}



template< typename SimplexTree, typename Filtration >
inline SimplexTree filtrationDionysusToGudhi(const Filtration & filtration) {

  std::map< typename Filtration::Simplex, unsigned,
      typename Filtration::Simplex::VertexComparison > simplex_map;
  unsigned size_of_simplex_map = 0;
  SimplexTree smplxTree;

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    std::vector< double > cmplxVec;
    double value;
    std::vector< double > boundaryVec;
    filtrationDionysusOne(c, simplex_map, 0, cmplxVec, value, boundaryVec);

    smplxTree.insert_simplex(cmplxVec, value);

    simplex_map.insert(typename std::map< typename Filtration::Simplex,
        unsigned >::value_type(c, size_of_simplex_map++));
  }

  return smplxTree;
}



template< typename Column, typename Dimension, typename Filtration,
          typename VectorList, typename RealVector, typename Boundary >
inline void filtrationDionysusToPhat(
    const Filtration & filtration, VectorList & cmplx, RealVector & values,
    Boundary & boundary_matrix) {

  const unsigned nFltr = filtration.size();
  std::map< typename Filtration::Simplex, typename Column::value_type,
    typename Filtration::Simplex::VertexComparison > simplex_map;
  typename Column::value_type size_of_simplex_map = 0;

  cmplx = VectorList(nFltr);
  values = RealVector(nFltr);
  boundary_matrix.set_num_cols(nFltr);
  typename VectorList::iterator iCmplx = cmplx.begin();
  typename RealVector::iterator iValue = values.begin();

  for (typename Filtration::Index it = filtration.begin();
       it != filtration.end(); ++it, ++iCmplx, ++iValue) {
    const typename Filtration::Simplex & c = filtration.simplex(it);

    Column cmplxVec;
    Column boundary_indices;
    filtrationDionysusOne(
        c, simplex_map, 0, cmplxVec, *iValue, boundary_indices);
    *iCmplx = cmplxVec;

    std::sort(boundary_indices.begin(), boundary_indices.end());
    boundary_matrix.set_col(size_of_simplex_map, boundary_indices);
    Dimension dim_of_column = c.dimension();
    boundary_matrix.set_dim(size_of_simplex_map, dim_of_column);

    simplex_map.insert(typename std::map< typename Filtration::Simplex,
        typename Column::value_type >::value_type(c, size_of_simplex_map++));
  }
}



# endif // __TYPECASTUTILS_H__
