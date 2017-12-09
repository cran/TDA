// for kernel density
#include <tdautils/kernelUtils.h>
#include <string>



// KDE function on a Grid
template< typename RealVector, typename RealMatrix, typename Print >
inline RealVector kde(
    const RealMatrix  & X,
    const RealMatrix  & Grid,
    const unsigned      nSample,
    const unsigned      dim,
    const unsigned      nGrid,
    const double        h,
    const std::string & kertype,
    const RealVector  & weight,
    const bool          printProgress,
    const Print       & print
) {
  const double pi = 3.141592653589793;
  const double h_to_dim = pow(h, (int)dim);
//  const double den = pow(h, (int)dim) * pow(2 * pi, dim / 2.0);
  RealVector kdeValue;
  int counter = 0, percentageFloor = 0;
  int totalCount = nGrid;

  if (printProgress) {
    printProgressFrame(print);
  }

  double(*kernel)(const double);

  if (dim <= 1 || kertype.length() > 12) {
    if (kertype[0] == 'G' || kertype[0] == 'g') {
      kernel = gaussianSquare;
    }
    if (kertype[0] == 'E' || kertype[0] == 'e') {
      kernel = epanechnikovSquare;
    }
    kdeValue = computeKernel< RealVector >(
        X, Grid, nSample, dim, nGrid, h * h, kernel, weight, printProgress,
        print, counter, totalCount, percentageFloor);
  }
  else {
    if (kertype[0] == 'G' || kertype[0] == 'g') {
      kernel = gaussian;
    }
    if (kertype[0] == 'E' || kertype[0] == 'e') {
      kernel = epanechnikov;
    }
    kdeValue = computeGaussOuter< RealVector >(
        X, Grid, nSample, dim, nGrid, h, kernel, weight, printProgress,
        print, counter, totalCount, percentageFloor);
  }

  for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
    kdeValue[gridIdx] /= h_to_dim;
  }

  if (printProgress) {
    print("\n");
  }

  return kdeValue;
}



// kernel Dist function on a Grid
template< typename RealVector, typename RealMatrix, typename Print >
inline RealVector kdeDist(
    const RealMatrix & X,
    const RealMatrix & Grid,
    const unsigned     nSample,
    const unsigned     dim,
    const unsigned     nGrid,
    const double       h,
    const RealVector & weight,
    const bool         printProgress,
    const Print      & print
) {

  // first = sum K_h(X_i, X_j), second = K_h(x, x), third = sum K_h(x, X_i)
  std::vector< double > firstValue;
  const double second = 1.0;
  std::vector< double > thirdValue;
  double firstmean;
  RealVector kdeDistValue(nGrid);
  int counter = 0, percentageFloor = 0;
  int totalCount = nSample + nGrid;

  if (printProgress) {
    printProgressFrame(print);
  }

  firstValue = computeKernel< std::vector< double > >(
      X, X, X.nrow(), X.ncol(), X.nrow(), h * h, gaussianSquare, weight,
      printProgress, print, counter, totalCount, percentageFloor);

  if (dim <= 1) {
    thirdValue = computeKernel< std::vector< double > >(
        X, Grid, X.nrow(), Grid.ncol(), Grid.nrow(), h * h, gaussianSquare,
        weight, printProgress, print, counter, totalCount, percentageFloor);
  }
  else {
    thirdValue = computeGaussOuter< std::vector< double > >(
        X, Grid, X.nrow(), Grid.ncol(), Grid.nrow(), h, gaussian, weight,
        printProgress, print, counter, totalCount, percentageFloor);
  }

  if (weight.size() == 1) {
    firstmean = std::accumulate(firstValue.begin(), firstValue.end(), 0.0) / nSample;
  }
  else {
    firstmean = std::inner_product(
      firstValue.begin(), firstValue.end(), weight.begin(), 0.0) /
      std::accumulate(weight.begin(), weight.end(), 0.0);
  }

  for (unsigned gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
    kdeDistValue[gridIdx] = std::sqrt(firstmean + second - 2 * thirdValue[gridIdx]);
  }

  if (printProgress) {
    print("\n");
  }

  return kdeDistValue;
}

