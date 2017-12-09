// distance to measure function on a Grid
template< typename RealVector, typename RealMatrix >
RealVector dtm(
    const RealMatrix & knnDistance,
    const unsigned     nGrid,
    const double       weightBound,
    const double       r
) {

  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
  Rcpp::NumericVector dtmValue(nGrid, 0.0);
  unsigned weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
        ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        dtmValue[gridIdx] += distanceTemp * distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
        (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
        ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        dtmValue[gridIdx] += distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp *
        (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
        ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        dtmValue[gridIdx] += std::pow(distanceTemp, r);
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
        (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}



// distance to measure function on a Grid, with weight
template< typename RealVector, typename RealMatrix >
RealVector dtmWeight(
    const RealMatrix & knnDistance,
    const unsigned     nGrid,
    const double       weightBound,
    const double       r,
    const RealMatrix & knnIndex,
    const RealVector & weight
) {

  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
  Rcpp::NumericVector dtmValue(nGrid, 0.0);
  double weightTemp, weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
           ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        weightTemp = weight[knnIndex[gridIdx + kIdx * nGrid] - 1];
        dtmValue[gridIdx] += distanceTemp * distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
        (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
           ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        weightTemp = weight[knnIndex[gridIdx + kIdx * nGrid] - 1];
        dtmValue[gridIdx] += distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * (weightBound - weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < nGrid; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
           ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * nGrid];
        weightTemp = weight[knnIndex[gridIdx + kIdx * nGrid] - 1];
        dtmValue[gridIdx] += std::pow(distanceTemp, r) * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
        (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}
