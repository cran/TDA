#include <cmath>
#include <numeric>
#include <map>

// read element of matrix
double ReadMat(double*XX, int *pNN, int *pDD, int i, int d){
	double out=0.0;
	out=XX[(d-1)*(pNN[0])+i-1];
	return out;
}


// write element of matrix
void WriteMat(double*XX, int *pNN, int *pDD, int i, int d, double input){
	XX[(d-1)*(pNN[0])+i-1]=input;
}


// print frame of progress
template <typename Print>
inline void printProgressFrame(Print print) {

	print("0   10   20   30   40   50   60   70   80   90   100");
	print("\n");
	print("|----|----|----|----|----|----|----|----|----|----|\n");
	print("*");
}


// print progress amount
template <typename Print>
inline void printProgressAmount(Print print, int& counter, const int totalCount, int& percentageFloor) {

	int progressAmount = std::floor((100 * (++counter) / totalCount - percentageFloor) / 2);
	if (progressAmount > 0) {
		for (int progressIdx = 1; progressIdx <= progressAmount; ++progressIdx) {
			print("*");
			percentageFloor += 2;
		}
	}
}


// get row of matrix
template <typename RealVector, typename RealMatrix>
inline RealVector matrixRow(
    const RealMatrix & X, const unsigned nX, const unsigned dim,
    const unsigned rowIdx) {
	RealVector rowVector(dim);
	for (unsigned colIdx = 0; colIdx < dim; ++colIdx) {
		rowVector[colIdx] = X[rowIdx + colIdx * nX];
	}
	return rowVector;
}



const double pi = 3.141592653589793;
const double one_over_sqrt_2_pi = 1.0 / sqrt(2.0 * pi);



// Gaussian Kernel
inline double gaussian(double x) {
  return one_over_sqrt_2_pi * exp(-x * x / 2);
}



// Gaussian Kernel
inline double gaussianSquare(double xSquare) {
  return one_over_sqrt_2_pi * exp(-xSquare / 2);
}



// Epanechnikov Kernel
inline double epanechnikov(double x) {
  double xSquare = x * x;
  if (xSquare < 1.0) {
    return 0.75 * (1.0 - xSquare);
  }
  else {
    return 0.0;
  }
}



// Epanechnikov Kernel
inline double epanechnikovSquare(double xSquare) {
  if (xSquare < 1.0) {
    return 0.75 * (1.0 - xSquare);
  }
  else {
    return 0.0;
  }
}



// oneKernel
template <typename RealVector1, typename RealVector2, typename RealMatrix>
inline double oneKernel(
    const RealVector1 & point, const RealMatrix & X, const unsigned nSample,
    const double hSquare, double (* kernelSquare)(const double),
    const RealVector2 & weight) {

	const unsigned dim = point.size();
	double sum, tmp;
	double oneKernelValue = 0.0;

	if (weight.size() == 1) {
		for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * nSample];
				sum += tmp * tmp;
			}
			oneKernelValue += kernelSquare(sum / hSquare);
		}
		return (oneKernelValue / nSample);

	}
	else {
		for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * nSample];
				sum += tmp * tmp;
			}
			oneKernelValue += kernelSquare(sum / hSquare) * weight[sampleIdx];
		}
		return (oneKernelValue / std::accumulate(weight.begin(), weight.end(), 0.0));
	}
}


// computes kernel with directly calculating ||X_i - Y_j||^2
template<typename RealVector1, typename RealMatrix1, typename RealMatrix2, typename RealVector2, typename Print>
inline RealVector1 computeKernel(
    const RealMatrix1 & X, const RealMatrix2 & Y, const unsigned nX,
    const unsigned dim, const unsigned nY, const double hSquare,
    double(*kernelSquare)(const double), const RealVector2 & weight,
    const bool printProgress, Print print, int & counter, int & totalCount,
    int & percentageFloor) {

	RealVector1 kernelValue(nY);

	if (printProgress) {

		for (unsigned yIdx = 0; yIdx < nY; ++yIdx) {
			kernelValue[yIdx] = oneKernel(
				  matrixRow< std::vector< double > >(Y, nY, dim, yIdx), X, nX, hSquare,
          kernelSquare, weight);

			// printProgress
			printProgressAmount(print, counter, totalCount, percentageFloor);
		}
	}
	else { //no printProgress
		for (unsigned yIdx = 0; yIdx < nY; ++yIdx) {
			kernelValue[yIdx] = oneKernel(
				  matrixRow< std::vector< double > >(Y, nY, dim, yIdx), X, nX, hSquare,
          kernelSquare, weight);
		}
	}

	return kernelValue;
}


// marginalize Grid
// GridMargin is a vector of marginal values of Grid
// marginIndex stores for each grid points and dimension
// correspond to which value in GridMargin
template <typename RealMatrix, typename RealVector, typename NonnegativeMatrix>
inline void marginalizeGrid(
    const RealMatrix & Grid, const unsigned dim, const unsigned nGrid,
    RealVector & GridMargin, NonnegativeMatrix & marginIndex) {

	marginIndex = NonnegativeMatrix(dim * nGrid);
	std::map< double, unsigned > GridMarginMap;
	for (unsigned dimIdx = 0, mgnIdx = 0; dimIdx < dim; ++dimIdx) {
		GridMarginMap.clear();
		for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
			if (GridMarginMap.find(Grid[gridIdx + dimIdx * nGrid]) ==
					GridMarginMap.end()) {
				GridMarginMap.insert(std::pair< double, unsigned >(
						Grid[gridIdx + dimIdx * nGrid], 0));
			}
		}
		GridMargin.resize(GridMargin.size() + GridMarginMap.size());
		for (std::map< double, unsigned >::iterator marginMapItr = GridMarginMap.begin();
		marginMapItr != GridMarginMap.end(); ++marginMapItr, ++mgnIdx) {
			GridMargin[mgnIdx] = marginMapItr->first;
			marginMapItr->second = mgnIdx;
		}
		for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
			marginIndex[dimIdx + gridIdx * dim] =
				GridMarginMap.find(Grid[gridIdx + dimIdx * nGrid])->second;
		}
	}
}


// gaussOuter
// returns matrix, with each element as O_ij = e(-1/(2h^2) (X_i - Y_j)^2) 
template< typename RealMatrix, typename RealVector1, typename RealVector2,
          typename Print >
inline RealMatrix GaussOuter(
    const RealVector1 & X, const RealVector2 & Y, const unsigned yNum,
    const double h, double(*kernel)(const double), bool printProgress,
    Print print, int & counter, const int totalCount, int & percentageFloor) {

	const unsigned xNum = X.size();
	//const unsigned yNum = Y.length();
	RealMatrix outerProd(xNum * yNum);
	if (printProgress) {
		for (unsigned xIdx = 0; xIdx < xNum; ++xIdx) {
			for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
        outerProd[yIdx + xIdx * yNum] = kernel((X[xIdx] - Y[yIdx]) / h);
					//exp(-(X[xIdx] - Y[yIdx]) * (X[xIdx] - Y[yIdx]) / (2 * h * h));
			}
			printProgressAmount(print, counter, totalCount, percentageFloor);
		}
	}
	else {
		for (unsigned xIdx = 0; xIdx < xNum; ++xIdx) {
			for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
				outerProd[yIdx + xIdx * yNum] = kernel((X[xIdx] - Y[yIdx]) / h);
//					exp(-(X[xIdx] - Y[yIdx]) * (X[xIdx] - Y[yIdx]) / (2 * h * h));
			}
		}
	}
	return outerProd;
}


// productCross
// cross product over different of matrices
template< typename RealVector1, typename RealMatrix,
          typename NonnegativeMatrix, typename RealVector2, typename Print >
inline RealVector1 productCross(const RealMatrix & outerValue,
		const NonnegativeMatrix & marginIndex, const RealVector2 & weight,
		const unsigned nSample, const unsigned dim, const unsigned nGrid,
		bool printProgress, Print print, int & counter, const int totalCount,
		int & percentageFloor) {

	RealVector1 gridValue(nGrid);
	if (printProgress) {
		if (weight.size() == 1) {
			for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * nSample +
							marginIndex[dimIdx + gridIdx * dim] * dim * nSample];
					}
					gridValue[gridIdx] += tmp;
				}
				gridValue[gridIdx] /= nSample;
				printProgressAmount(print, counter, totalCount, percentageFloor);
			}
		}
		else {
			const unsigned weightSum = std::accumulate(weight.begin(), weight.end(), 0.0);
			for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * nSample +
							marginIndex[dimIdx + gridIdx * dim] * nSample * dim];
					}
					gridValue[gridIdx] += tmp * weight[sampleIdx];
				}
				gridValue[gridIdx] /= weightSum;
				printProgressAmount(print, counter, totalCount, percentageFloor);
			}
		}
	}
	else {
		if (weight.size() == 1) {
			for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * nSample +
							marginIndex[dimIdx + gridIdx * dim] * dim * nSample];
					}
					gridValue[gridIdx] += tmp;
				}
				gridValue[gridIdx] /= nSample;
			}
		}
		else {
			const unsigned weightSum = std::accumulate(weight.begin(), weight.end(), 0.0);
			for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < nSample; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * nSample +
							marginIndex[dimIdx + gridIdx * dim] * nSample * dim];
					}
					gridValue[gridIdx] += tmp * weight[sampleIdx];
				}
				gridValue[gridIdx] /= weightSum;
			}
		}
	}
	return gridValue;
}


// computes Gaussian kernel
// with factorizing as ||X_i - Y_j||^2 = sum (X_ik - Y_jk)^2
// and using outer products to aggregate them
template< typename RealVector1, typename RealMatrix, typename RealVector2,
          typename Print >
inline RealVector1 computeGaussOuter(
    const RealMatrix & X, const RealMatrix & Grid, const unsigned nSample,
    const unsigned dim, const unsigned nGrid, const double h,
    double(*kernel)(const double), const RealVector2 & weight,
    const bool printProgress, Print print, int & counter, int & totalCount,
    int & percentageFloor) {
	
	std::vector< double > GridMargin, gaussOuter;
	std::vector< unsigned > marginIndex;
	RealVector1 gaussValue(nGrid);

	marginalizeGrid(Grid, dim, nGrid, GridMargin, marginIndex);

	totalCount += GridMargin.size();

	gaussOuter = GaussOuter< std::vector< double > >(
			GridMargin, X, nSample * dim, h, kernel, printProgress, print,
      counter, totalCount, percentageFloor);

	gaussValue = productCross< RealVector1 >(
			gaussOuter, marginIndex, weight, nSample, dim, nGrid,
			printProgress, print, counter, totalCount, percentageFloor);

	return gaussValue;
}
