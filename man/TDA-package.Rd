\name{TDA-package}

\alias{TDA-package}
\alias{TDA}

\docType{package}

\title{
  Statistical Tools for Topological Data Analysis
}

\description{
  Tools for Topological Data Analysis. The package focuses on statistical analysis of persistent homology and density clustering. For that, this package provides an R interface for the efficient algorithms of the C++ libraries \href{ https://project.inria.fr/gudhi/software/ }{GUDHI}, \href{ https://www.mrzv.org/software/dionysus/ }{Dionysus} and \href{ https://bitbucket.org/phat-code/phat/ }{PHAT}. This package also implements methods from Fasy et al. (2014) and Chazal et al. (2015) for analyzing the statistical significance of persistent homology features.
}

\details{
\tabular{ll}{
Package: \tab TDA\cr
Version: \tab 1.9.4\cr
Date: \tab 2025-02-01\cr
License: \tab GPL-3 \cr
}
}

\author{
  Brittany Terese Fasy, Jisu Kim, Fabrizio Lecci, Clement Maria, David L. Millman, and Vincent Rouvreau

  Maintainer: Jisu Kim <jkim82133@snu.ac.kr>
}

\references{
  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.

  Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2014). "Confidence Sets for Persistence Diagrams." Annals of Statistics. (arXiv:1303.7117)

  Chazal F, Fasy BT, Lecci F, Rinaldo A, Wasserman L (2015). "Stochastic Convergence of Persistence Landscapes and Silhouettes." Journal of Computational Geometry. (arXiv:1312.0308)

  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2015a). "Subsampling Methods for Persistent Homology." Proceedings of the 32nd International Conference on Machine Learning (ICML). (arXiv:1406.1901)

  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2017). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Journal of Machine Learning Research. (arXiv:1412.7197)

  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology." \url{ https://www.mrzv.org/software/dionysus/ }.

  Bauer U, Kerber M, Reininghaus J (2012). "PHAT, a software library for persistent homology". \url{ https://bitbucket.org/phat-code/phat/ }.
}

\keyword{package}
