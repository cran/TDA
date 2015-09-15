#include <cmath>

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
inline RealVector matrixRow(const RealMatrix& X, const unsigned rowIdx) {
	const unsigned colNum = X.ncol();
	const unsigned rowNum = X.nrow();
	RealVector rowVector(colNum);
	for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
		rowVector[colIdx] = X[rowIdx + colIdx * rowNum];
	}
	return rowVector;
}


// oneKernel
template <typename RealVector1, typename RealVector2, typename RealMatrix>
inline double oneKernel(const RealVector1& point, const RealMatrix& X, const double h, const RealVector2& weight) {

	const unsigned dimension = X.ncol();
	const unsigned sampleNum = X.nrow();
	double sum, tmp;
	double oneKernelValue = 0.0;

	if (weight.size() == 1) {
		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * sampleNum];
				sum += tmp * tmp;
			}
			oneKernelValue += exp(-sum / (2 * h * h));
		}
		return (oneKernelValue / sampleNum);

	}
	else {
		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * sampleNum];
				sum += tmp * tmp;
			}
			oneKernelValue += exp(-sum / (2 * h * h)) * weight[sampleIdx];
		}
		return (oneKernelValue / std::accumulate(weight.begin(), weight.end(), 0.0));
	}
}


// computes kernel with directly calculating ||X_i - Y_j||^2
template<typename RealVector1, typename RealMatrix1, typename RealMatrix2, typename RealVector2, typename Print>
inline RealVector1 computeKernel(const RealMatrix1& X,
	const RealMatrix2& Y, const double h, const RealVector2& weight,
	const bool printProgress, Print print, int& counter, int& totalCount,
	int& percentageFloor) {

	const unsigned yNum = Y.nrow();
	RealVector1 kernelValue(yNum);

	if (printProgress) {
		for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
			kernelValue[yIdx] = oneKernel(
				matrixRow< std::vector< double > >(Y, yIdx), X, h, weight);

			// printProgress
			printProgressAmount(print, counter, totalCount, percentageFloor);
		}
	}
	else { //no printProgress
		for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
			kernelValue[yIdx] = oneKernel(
				matrixRow< std::vector< double > >(Y, yIdx), X, h, weight);
		}
	}

	return kernelValue;
}


// marginalize Grid
// GridMargin is a vector of marginal values of Grid
// marginIndex stores for each grid points and dimension
// correspond to which value in GridMargin
template <typename RealMatrix, typename RealVector, typename NonnegativeMatrix>
inline void marginalizeGrid(const RealMatrix& Grid, RealVector& GridMargin,
		NonnegativeMatrix& marginIndex) {

	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	marginIndex = NonnegativeMatrix(dimension * gridNum);
	std::map< double, unsigned > GridMarginMap;
	for (unsigned dimIdx = 0, mgnIdx = 0; dimIdx < dimension; ++dimIdx) {
		GridMarginMap.clear();
		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			if (GridMarginMap.find(Grid[gridIdx + dimIdx * gridNum]) ==
					GridMarginMap.end()) {
				GridMarginMap.insert(std::pair< double, unsigned >(
						Grid[gridIdx + dimIdx * gridNum], 0));
			}
		}
		GridMargin.resize(GridMargin.size() + GridMarginMap.size());
		for (std::map< double, unsigned >::iterator marginMapItr = GridMarginMap.begin();
		marginMapItr != GridMarginMap.end(); ++marginMapItr, ++mgnIdx) {
			GridMargin[mgnIdx] = marginMapItr->first;
			marginMapItr->second = mgnIdx;
		}
		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			marginIndex[dimIdx + gridIdx * dimension] =
				GridMarginMap.find(Grid[gridIdx + dimIdx * gridNum])->second;
		}
	}
}


// gaussOuter
// returns matrix, with each element as O_ij = e(-1/(2h^2) (X_i - Y_j)^2) 
template <typename RealMatrix, typename RealVector1, typename RealVector2, typename Print>
inline RealMatrix GaussOuter(const RealVector1& X, const RealVector2& Y,
		const double h, bool printProgress, Print print, int& counter,
		const int totalCount, int& percentageFloor) {

	const unsigned xNum = X.size();
	const unsigned yNum = Y.length();
	RealMatrix outerProd(xNum * yNum);
	if (printProgress) {
		for (unsigned xIdx = 0; xIdx < xNum; ++xIdx) {
			for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
				outerProd[yIdx + xIdx * yNum] =
					exp(-(X[xIdx] - Y[yIdx]) * (X[xIdx] - Y[yIdx]) / (2 * h * h));
			}
			printProgressAmount(print, counter, totalCount, percentageFloor);
		}
	}
	else {
		for (unsigned xIdx = 0; xIdx < xNum; ++xIdx) {
			for (unsigned yIdx = 0; yIdx < yNum; ++yIdx) {
				outerProd[yIdx + xIdx * yNum] =
					exp(-(X[xIdx] - Y[yIdx]) * (X[xIdx] - Y[yIdx]) / (2 * h * h));
			}
		}
	}
	return outerProd;
}


// productCross
// cross product over different of matrices
template <typename RealVector1, typename RealMatrix, typename NonnegativeMatrix, typename RealVector2, typename Print>
inline RealVector1 productCross(const RealMatrix& outerValue,
		const NonnegativeMatrix& marginIndex, const RealVector2& weight,
		const unsigned sampleNum, const unsigned dimension, const unsigned gridNum,
		bool printProgress, Print print, int& counter, const int totalCount,
		int& percentageFloor) {

	RealVector1 gridValue(gridNum);
	if (printProgress) {
		if (weight.size() == 1) {
			for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * sampleNum +
							marginIndex[dimIdx + gridIdx * dimension] * dimension * sampleNum];
					}
					gridValue[gridIdx] += tmp;
				}
				gridValue[gridIdx] /= sampleNum;
				printProgressAmount(print, counter, totalCount, percentageFloor);
			}
		}
		else {
			const unsigned weightSum = std::accumulate(weight.begin(), weight.end(), 0.0);
			for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * sampleNum +
							marginIndex[dimIdx + gridIdx * dimension] * sampleNum * dimension];
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
			for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * sampleNum +
							marginIndex[dimIdx + gridIdx * dimension] * dimension * sampleNum];
					}
					gridValue[gridIdx] += tmp;
				}
				gridValue[gridIdx] /= sampleNum;
			}
		}
		else {
			const unsigned weightSum = std::accumulate(weight.begin(), weight.end(), 0.0);
			for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
				gridValue[gridIdx] = 0.0;
				for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
					double tmp = 1.0;
					for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
						tmp *= outerValue[sampleIdx + dimIdx * sampleNum +
							marginIndex[dimIdx + gridIdx * dimension] * sampleNum * dimension];
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
template<typename RealVector1, typename RealMatrix1, typename RealMatrix2, typename RealVector2, typename Print>
inline RealVector1 computeGaussOuter(const RealMatrix1& X,
		const RealMatrix2& Grid, const double h, const RealVector2& weight,
		const bool printProgress, Print print, int& counter, int& totalCount,
		int& percentageFloor) {
	
	const unsigned sampleNum = X.nrow();
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	std::vector< double > GridMargin, gaussOuter;
	std::vector< unsigned > marginIndex;
	RealVector1 gaussValue(gridNum);

	marginalizeGrid(Grid, GridMargin, marginIndex);

	totalCount += GridMargin.size();

	gaussOuter = GaussOuter< std::vector< double > >(
			GridMargin, X, h, printProgress, print, counter, totalCount,
			percentageFloor);

	gaussValue = productCross< RealVector1 >(
			gaussOuter, marginIndex, weight, sampleNum, dimension, gridNum,
			printProgress, print, counter, totalCount, percentageFloor);

	return gaussValue;
}
