#include <topology/rips.h>
#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>

#include <geometry/l2distance.h>
#include <geometry/Arbitdistance.h>
#include <geometry/distances.h>
#include <utilities/containers.h>           // for BackInsertFunctor
#include <utilities/timer.h>

#include <vector>


typedef         PairwiseDistances<PointContainer, ArbitDistance>        PairDistancesA;
typedef         PairDistancesA::DistanceType                             DistanceTypeA;
typedef         PairDistancesA::IndexType                                VertexRA;
typedef         Rips< PairDistancesA, Simplex< VertexRA, double > >      GeneratorA;
typedef         GeneratorA::Simplex                                      SmplxRA;
typedef         Filtration<SmplxRA>                                       FltrRA;
typedef         StaticPersistence<>                                     PersistenceR;
//typedef         DynamicPersistenceChains<>                              PersistenceR;
typedef         PersistenceDiagram<>                                    PDgmR;
