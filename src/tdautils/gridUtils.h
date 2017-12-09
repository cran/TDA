#ifndef __GRIDUTILS_H__
#define __GRIDUTILS_H__

#include <utilities/log.h>

#include <topology/simplex.h>
#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>
#include <utilities/indirect.h>

#include <vector>
#include <map>
#include <numeric>


#if 1
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif



typedef         unsigned                                            Vertex;
typedef         Simplex<Vertex, double>                             Smplx;
typedef         Smplx::VertexContainer				    VertexCont;
// typedef         std::vector<Vertex>                                 VertexVector;
typedef         Filtration<Smplx>                                   Fltr;
typedef         StaticPersistence<>                                 Persistence;
typedef         PersistenceDiagram<>                                PDgm;
typedef         OffsetBeginMap<Persistence, Fltr, 
                               Persistence::iterator, 
                               Fltr::Index>                         PersistenceFiltrationMap;
typedef         OffsetBeginMap<Fltr, Persistence,
                               Fltr::Index, 
                               Persistence::iterator>               FiltrationPersistenceMap;



// add a single edge to the filtration
template< typename VectorList >
void addEdge(int vert01, int vert02, VectorList & cmplx) {
     typename VectorList::value_type vertices(2);
     vertices[0] = vert01;
     vertices[1] = vert02;
     cmplx.push_back(vertices);
} // end function to add a single edge



// add a single triangle to the filtration
template< typename VectorList >
void addTri(int vert01, int vert02, int vert03, VectorList & cmplx) {
     typename VectorList::value_type vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     cmplx.push_back(vertices);
} // end function to add a single triangle



// add a single tet to the filtration
template< typename VectorList >
void addTet(int vert01, int vert02, int vert03, int vert04, VectorList & cmplx) {
     typename VectorList::value_type vertices(4);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     vertices[3] = vert04;
     cmplx.push_back(vertices);
} // end function to add a single tet



template< typename VectorList >
void addAllEdges(
    const int ncols, const int nrows, int i, int j, int k,
    VectorList & cmplx) {     
     int curidx = i + ncols*j + ncols*nrows*k;

     // ... add edge (i-1,j,k) <--> (i,j,k)
     if (i > 0)
     {
        addEdge(curidx, curidx -1, cmplx);
     }
     
     // ... add edge (i,j-1,k) <--> (i,j,k)
     if (j > 0)
     {
        addEdge(curidx, curidx - ncols, cmplx);
     }
   
     // ... add edge (i,j,k-1) <--> (i,j,k)
     if (k > 0)
     {
        addEdge(curidx, curidx - nrows*ncols, cmplx);
     }


     // TODO: add the rest of the code for creating edges to here
     //       from fcn simplicesFromGrid 

     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX 
        if (i > 0 && j > 0) // top
        { 
           addEdge(curidx, curidx - ncols -1, cmplx);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(curidx, curidx - nrows*ncols -1, cmplx);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(curidx, curidx - nrows*ncols - ncols, cmplx);
        }
     }
     else
     {
     	// ... ODD BOX
        if (i > 0 && j > 0) // top
        {
           addEdge(curidx - 1, curidx - ncols, cmplx);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(curidx - 1, curidx - nrows*ncols, cmplx);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(curidx - ncols, curidx - nrows*ncols, cmplx);
        }
     }


     return;
} // end function addEdges



template< typename VectorList >
void addEvenTets(
    const int ncols, const int nrows, int i, int j, int k,
    VectorList & cmplx) {

     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // top vertex (i, j-1, k)
     addTet(curidx, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx - ncols, cmplx);
     
     // top vertex (i-1, j, k)
     addTet(curidx, curidx - 1, curidx - nrows*ncols - 1, curidx -1 - ncols, cmplx);

     // top vertex (i, j, k-1)
     addTet(curidx, curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - nrows*ncols, cmplx);

     // top vertex (i-1, j-1, k-1)
     addTet(curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - 1 - ncols, curidx - 1 - ncols - nrows*ncols, cmplx);
     
     return;
} // end fcn to add four EVEN tets



template< typename VectorList >
void addOddTets(
    const int ncols, const int nrows, int i, int j, int k,
    VectorList & cmplx) {
     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
      
     typename VectorList::value_type vertices(4);
     vertices[0] = curidx;  vertices[3] = curidx;
     vertices[1] = -1; vertices[2] = -1; 
    
     int v1, v2, v3, v4;
     double value, value2;  // max of value and value 2 is the fcn value. 

     // top vertex (i, j, k)
     v1 = curidx -1;   vertices[0] = v1;
     v2 = curidx - ncols;  vertices[1] = v2;  
     
     v3 = curidx - nrows*ncols;  vertices[2] = v3;
     v4 = curidx; vertices[3] = v4;
     
     cmplx.push_back(vertices);
     
     // top vertex (i-1, j-1, k)
     v3 = curidx - 1 - ncols - nrows*ncols;  vertices[2] = v3;
     v4 = curidx - 1 -ncols; vertices[3] = v4;
     
     cmplx.push_back(vertices);

     // top vertex (i, j-1, k-1)
     v1 = curidx - nrows*ncols;   vertices[0] = v1;
     v2 = curidx - 1 - ncols - nrows*ncols;  vertices[1] = v2;  
     
     v3 = curidx - ncols;  vertices[2] = v3;
     v4 = curidx -ncols - nrows*ncols; vertices[3] = v4;
     
     cmplx.push_back(vertices);

     // top vertex (i-1, j, k-1)
     v3 = curidx - 1;  vertices[2] = v3;
     v4 = curidx -1 - nrows * ncols; vertices[3] = v4;
    
     cmplx.push_back(vertices);

} // end fcn addEvenTets



template< typename VectorList >
void addAllTriangles(
    const int ncols, const int nrows, int i, int j, int k,
    VectorList & cmplx) {

     int curidx = i + ncols*j + ncols*nrows*k;
     
     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX
        if (i > 0 && j > 0) // top
        {
     
           addTri(curidx, curidx - ncols - 1, curidx - ncols, cmplx);
           addTri(curidx, curidx -1, curidx - ncols -1, cmplx); 
        }
        if (i > 0 && k > 0) // back
        {
           addTri(curidx, curidx - nrows*ncols -1, curidx - 1, cmplx);
           addTri(curidx, curidx - nrows*ncols, curidx - nrows*ncols -1, cmplx);
        }
        
        if (j > 0 && k > 0) // right
        {
           addTri(curidx, curidx - nrows*ncols - ncols, curidx - nrows*ncols, cmplx);
           addTri(curidx, curidx - ncols, curidx - nrows*ncols - ncols, cmplx);
           
           if (i > 0) // middle
           {
              addTri(curidx, curidx - ncols -1, curidx - ncols - nrows*ncols, cmplx);
              addTri(curidx, curidx -1 -nrows*ncols, curidx - ncols -1, cmplx);
              addTri(curidx -1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx, cmplx);
              addTri(curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols, cmplx);
           }
        }
     } // end if for even case 
     else {
        // ... ODD CASE
        if (i > 0 && j > 0) // top
        {
           addTri(curidx -1, curidx - ncols, curidx, cmplx);
           addTri(curidx -1, curidx - ncols -1, curidx - ncols, cmplx);
        }

        if (i > 0 && k > 0) // back
        {
           addTri(curidx -1, curidx - nrows*ncols, curidx - nrows*ncols - 1, cmplx);
           addTri(curidx -1, curidx, curidx - nrows*ncols, cmplx);
        }  
        
        if (j > 0 && k > 0) // right
        {
           addTri(curidx - ncols, curidx - nrows*ncols, curidx - ncols - nrows*ncols, cmplx);
           addTri(curidx - ncols, curidx, curidx - nrows*ncols, cmplx);
           
           if (i > 0) // middle
           { 
              addTri(curidx -1, curidx - ncols, curidx - nrows*ncols, cmplx);
              addTri(curidx -1, curidx - nrows*ncols - ncols - 1, curidx - ncols, cmplx);
              addTri(curidx - nrows*ncols, curidx - nrows*ncols - ncols -1, curidx - ncols, cmplx);
              addTri(curidx - nrows*ncols, curidx - 1, curidx - nrows*ncols - ncols -1, cmplx);
           }
        } // end for through j k positive
     } // end else through odd case.

    return;
} // end function addAllTriangles



template< typename VectorList >
void addAllTetrahedra(
    const int ncols, const int nrows, int i, int j, int k,
    VectorList & cmplx) {
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX
        if (i > 0 && j > 0 && k > 0) // middle
        {
            // ... add center tets 
			addTet(curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx, cmplx);
            // ... add remaining tets 
			addEvenTets(ncols, nrows, i, j, k, cmplx);
        }
     } // end if for even case 
     else {
        // ... ODD CASE
        if (i > 0 && j > 0 && k > 0) // middle
        {
			// ... add central tet
			addTet(curidx -1, curidx - ncols, curidx - nrows*ncols, curidx - nrows*ncols -ncols -1, cmplx);
			// ... add remaining tets
			addOddTets(ncols, nrows, i, j, k, cmplx);
        } // end for through j k positive
     } // end else through odd case.

    return;
} // end function addTriangles



template< typename DimensionVector, typename VectorList >
void simplicesFromGrid(
    const DimensionVector & gridDim, const int embedDim, VectorList & cmplx) {

	const unsigned gridProd = std::accumulate(
			gridDim.begin(), gridDim.end(), 1, std::multiplies< int >());
	int ncols, nrows;
	ncols = nrows = 1;
	int i = 0; // indexing the columns
	int j = 0; // indexing the rows
	int k = 0; // indexing the z dimension
	unsigned int curidx = 0; // curidx = i + ncols * j + nrows * ncols * k

	if (gridDim.size() > 0) {
		ncols = gridDim[0];
	}
	if (gridDim.size() > 1) {
		nrows = gridDim[1];
	}

  while(curidx < gridProd) {

    // .. add the vertex 
    typename VectorList::value_type vcont;
    vcont.push_back((Vertex)(curidx));
    cmplx.push_back(vcont);

    // If dimension of embedded space >= 1, add the edges:
		if (embedDim >= 1) {
      addAllEdges(ncols, nrows, i, j, k, cmplx);
		}

		// If dimension of embedded space >= 2, add the triangles:
		if (embedDim >= 2) {
			addAllTriangles(ncols, nrows, i, j, k, cmplx);
		}

    // If dimension of embedded space >= 3, add the tetrahedra:
		if (embedDim >= 3) {
			addAllTetrahedra(ncols, nrows, i, j, k, cmplx);
      addAllTetrahedra(ncols, nrows, i, j, k, cmplx);
		}

    ++i; // advance column
    // advance row
    if (i >= ncols) {
      i = 0;
      ++j;
    }
    // advance z value
    if (j >= nrows) {
      j = 0;
      ++k;
    }
    ++curidx; // advance curidx

  }
} // end simplicesFromGrid function



template <typename IntVector>
inline std::vector< unsigned char > isInternal(unsigned int argIdx, const IntVector& gridDim) {
    std::vector< unsigned char > resIsInt;
    resIsInt.reserve(gridDim.size());
    typename IntVector::const_iterator itrDim;
    for (itrDim = gridDim.begin(); itrDim != gridDim.end(); itrDim++)
    {
        resIsInt.push_back((unsigned char)(argIdx % (*itrDim) > 0));
        argIdx /= (*itrDim);
    }
    return resIsInt;
}



inline std::vector< std::vector< unsigned char > > verticesLessVertex(const std::vector< unsigned char > & argVtx, const bool argAlsoEqual) {
	std::vector< std::vector< unsigned char > > resCubeVertices;
    unsigned int idxVtx, vtxNum;
    std::vector< unsigned int > oneTwoVec;
	oneTwoVec.reserve(argVtx.size());
	std::vector< unsigned char >::const_iterator itrVtx;
	for (itrVtx = argVtx.begin(); itrVtx != argVtx.end(); ++itrVtx)
	{
		oneTwoVec.push_back(1+(unsigned int)(*itrVtx));
	}
	
    vtxNum = std::accumulate(oneTwoVec.begin(), oneTwoVec.end(), 1, std::multiplies< unsigned int >());
	if (!argAlsoEqual)
	{
		vtxNum -= 1;
	}
    resCubeVertices.reserve(vtxNum);
    for (idxVtx = 0; idxVtx < vtxNum; ++idxVtx)
    {
        resCubeVertices.push_back(isInternal(idxVtx, oneTwoVec));
    }
	return resCubeVertices;
}



std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > triangulateHypercube(const int argDimEmbed, const unsigned char embedDim) {
    std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > resTriedCube;
    resTriedCube.reserve(embedDim+1);

    // vertices of hypercube
	std::vector< unsigned char > rootVtx(argDimEmbed, 1);
    std::vector< std::vector< unsigned char > > cubeVertices = verticesLessVertex(rootVtx, true);

	std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > mapDirSmpxVec;
	std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVec;
	std::vector< std::vector< unsigned char > > dirSmpx;
	std::vector< std::vector< unsigned char > >::const_iterator itrVtx;

	// 0 dim
	mapDirSmpxVec.clear();
	for (itrVtx = cubeVertices.begin(); itrVtx != cubeVertices.end(); ++itrVtx)
	{
		dirSmpxVec.clear();
		dirSmpx.clear();
		dirSmpx.push_back(*itrVtx);
		dirSmpxVec.push_back(dirSmpx);
		mapDirSmpxVec[ *itrVtx ] = dirSmpxVec;
	}
	resTriedCube.push_back(mapDirSmpxVec);

	unsigned char idxDim;
	std::vector< std::vector< unsigned char > > vtxLessVtx;
	std::vector< std::vector< unsigned char > >::const_iterator itrLessVtx;
	std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVecPrev;
	std::vector< std::vector< std::vector< unsigned char > > >::iterator itrSmpxVec;

	for (idxDim = 1; idxDim <= embedDim; ++idxDim)
	{
		mapDirSmpxVec.clear();
		for (itrVtx = cubeVertices.begin(); itrVtx != cubeVertices.end(); ++itrVtx)
		{
			dirSmpxVec.clear();
			vtxLessVtx = verticesLessVertex(*itrVtx, false);
			for (itrLessVtx = vtxLessVtx.begin(); itrLessVtx != vtxLessVtx.end(); ++itrLessVtx)
			{
				dirSmpxVecPrev = resTriedCube.at(idxDim-1).at(*itrLessVtx);
				for (itrSmpxVec = dirSmpxVecPrev.begin(); itrSmpxVec != dirSmpxVecPrev.end(); ++itrSmpxVec)
				{
					itrSmpxVec->push_back(*itrVtx); 
				}
				dirSmpxVec.insert(dirSmpxVec.end(), dirSmpxVecPrev.begin(), dirSmpxVecPrev.end());
			}
			mapDirSmpxVec[ *itrVtx ] = dirSmpxVec;
		}
		resTriedCube.push_back(mapDirSmpxVec);
	}

    return resTriedCube;
}



template< typename DimensionVector, typename VectorList >
void addSimplices(
    const int argIdxCur, const DimensionVector & gridDim, const unsigned char argIdxDim,
    std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > > & argTriedCube,
    VectorList & cmplx) {
    std::vector< unsigned char > isInt = isInternal(argIdxCur, gridDim);
    std::vector< std::vector< std::vector< unsigned char > > > dirSmpxVec = (argTriedCube.at(argIdxDim)).at(isInt);
    std::vector< std::vector< std::vector< unsigned char > > >::const_iterator itrDirSmpxVec;
    std::vector< std::vector< unsigned char > >::const_iterator itrDirVtxVec;
	std::vector< unsigned char > diffVtx(gridDim.size());
	std::vector< unsigned int > gridAccNum(gridDim.size(),1);
	std::partial_sum(gridDim.begin(), gridDim.end()-1, gridAccNum.begin()+1, std::multiplies< unsigned int >());

    typename VectorList::value_type vtxVec(argIdxDim + 1);
    typename VectorList::value_type::iterator itrVtxVec;
    for (itrDirSmpxVec = dirSmpxVec.begin(); itrDirSmpxVec != dirSmpxVec.end(); ++itrDirSmpxVec)
    {
        for (itrDirVtxVec = itrDirSmpxVec->begin(), itrVtxVec = vtxVec.begin();
             itrDirVtxVec != itrDirSmpxVec->end(); 
             ++itrDirVtxVec, ++itrVtxVec)
        {
			std::transform(isInt.begin(), isInt.end(), itrDirVtxVec->begin(), diffVtx.begin(), std::minus< char >());

            (*itrVtxVec) = (argIdxCur - std::inner_product(gridAccNum.begin(), gridAccNum.end(), diffVtx.begin(), 0));
        }
        cmplx.push_back(vtxVec);
    }
}



template< typename DimensionVector, typename VectorList >
void simplicesFromGridBarycenter(
    const DimensionVector & gridDim, const unsigned char embedDim,
    VectorList & cmplx) {

	const unsigned gridProd = std::accumulate(
			gridDim.begin(), gridDim.end(), 1, std::multiplies< int >());
	unsigned int idxCur; unsigned char idxDim;

   std::vector< std::map< std::vector< unsigned char >, std::vector< std::vector< std::vector< unsigned char > > > > >  triedCube = triangulateHypercube(gridDim.size(), embedDim);
  
    for (idxCur = 0; idxCur < gridProd ; ++idxCur) {
        for (idxDim = 0; idxDim <= embedDim; ++idxDim) {
		    addSimplices(idxCur, gridDim, idxDim, triedCube, cmplx);
        }
    }
}



# endif // __GRIDUTILS_H__
