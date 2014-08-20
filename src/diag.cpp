#include <R.h>
#include <R_ext/Print.h>



//for kernel density
#include <kernelUtilities.h>


//for bottleneck distance
#include <utilities/types.h>
#include <string>
#include <sstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/map.hpp>



// For rips and ripsArbit
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


typedef         PairwiseDistances<PointContainer, L2Distance>           PairDistances;
typedef         PairwiseDistances<PointContainer, ArbitDistance>        PairDistancesA;
typedef         PairDistances::DistanceType                             DistanceType;
typedef         PairDistances::IndexType                                VertexR;
typedef         PairDistancesA::DistanceType                             DistanceTypeA;
typedef         PairDistancesA::IndexType                                VertexRA;


typedef         Rips<PairDistances>                                     Generator;
typedef         Rips<PairDistancesA>                                     GeneratorA;

typedef         Generator::Simplex                                      SmplxR;
typedef         GeneratorA::Simplex                                      SmplxRA;

typedef         Filtration<SmplxR>                                       FltrR;
typedef         Filtration<SmplxRA>                                       FltrRA;

typedef         StaticPersistence<>                                     PersistenceR;
//typedef         DynamicPersistenceChains<>                              PersistenceR;
typedef         PersistenceDiagram<>                                    PDgmR;



// For grid
#include "utilities/log.h"
#include "topology/simplex.h"
#include "utilities/indirect.h"
#include <map>
#include <iostream>

#if 1
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif

typedef         unsigned                                            Vertex;
typedef         Simplex<Vertex, double>                             Smplx;
typedef         Smplx::VertexContainer				    VertexCont;
typedef         std::vector<Vertex>                                 VertexVector;
typedef         Filtration<Smplx>                                   Fltr;
typedef         StaticPersistence<>                                 Persistence;
typedef         PersistenceDiagram<>                                PDgm;
typedef         OffsetBeginMap<Persistence, Fltr, 
                               Persistence::iterator, 
                               Fltr::Index>                         PersistenceFiltrationMap;
typedef         OffsetBeginMap<Fltr, Persistence,
                               Fltr::Index, 
                               Persistence::iterator>               FiltrationPersistenceMap;










//bottleneck
typedef PersistenceDiagram<>                    PDgmB;

void read_diagram(PDgmB& dgm, const std::string& filename)
{
    std::ifstream in(filename.c_str());
    std::string line;
    std::getline(in, line);
    while (in)
    {
        std::istringstream sin(line);
        double x,y;
        sin >> x >> y;
        dgm.push_back(PDgmB::Point(x,y));
        std::getline(in, line);
    }
}





// add a single edge to the filtration
void addEdge(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(fcnvalues.at(vert01), fcnvalues.at(vert02));
     filtr.push_back(Smplx(bg, bg + 2, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single edge

// add a single triangle to the filtration
void addTri(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02, int vert03)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(std::max(fcnvalues.at(vert01), fcnvalues.at(vert02)),
                    fcnvalues.at(vert03));
     filtr.push_back(Smplx(bg, bg + 3, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single triangle

// add a single tet to the filtration
void addTet(Fltr& filtr, const std::vector<double> fcnvalues, 
            int vert01, int vert02, int vert03, int vert04)
{
     VertexVector vertices(3);
     vertices[0] = vert01;
     vertices[1] = vert02;
     vertices[2] = vert03;
     vertices[3] = vert04;
     VertexVector::const_iterator bg = vertices.begin();

     double value = std::max(std::max(fcnvalues.at(vert01), fcnvalues.at(vert02)),
                    std::max(fcnvalues.at(vert03), fcnvalues.at(vert04)));
     filtr.push_back(Smplx(bg, bg + 4, value)); 
         // std::max(fcnvalues.at(vert03),fcnvalues.at(vert04))),
} // end function to add a single tet

void addAllEdges(Fltr& filtr, const std::vector<double> fcnvalues, 
              const int ncols, const int nrows, int i, int j, int k)
{     
     int curidx = i + ncols*j + ncols*nrows*k;

     // ... add edge (i-1,j,k) <--> (i,j,k)
     if (i > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx -1);
     }
     
     // ... add edge (i,j-1,k) <--> (i,j,k)
     if (j > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx - ncols);
     }
   
     // ... add edge (i,j,k-1) <--> (i,j,k)
     if (k > 0)
     {
        addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols);
     }


     // TODO: add the rest of the code for creating edges to here
     //       from fcn simplicesFromGrid 

     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX 
        if (i > 0 && j > 0) // top
        { 
           addEdge(filtr, fcnvalues, curidx, curidx - ncols -1);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols -1);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(filtr, fcnvalues, curidx, curidx - nrows*ncols - ncols);
        }
     }
     else
     {
     	// ... ODD BOX
        if (i > 0 && j > 0) // top
        {
           addEdge(filtr, fcnvalues, curidx - 1, curidx - ncols);
        }
        if (i > 0 && k > 0) // back
        {
           addEdge(filtr, fcnvalues, curidx - 1, curidx - nrows*ncols);
        }
        if (j > 0 && k > 0) // right
        {
           addEdge(filtr, fcnvalues, curidx - ncols, curidx - nrows*ncols);
        }
     }


     return;
} // end function addEdges

void addEvenTets(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // top vertex (i, j-1, k)
     addTet(filtr, fcnvalues, curidx, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx - ncols);
     
     // top vertex (i-1, j, k)
     addTet(filtr, fcnvalues, curidx, curidx - 1, curidx - nrows*ncols - 1, curidx -1 - ncols);

     // top vertex (i, j, k-1)
     addTet(filtr, fcnvalues, curidx, curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - nrows*ncols);

     // top vertex (i-1, j-1, k-1)
     addTet(filtr, fcnvalues, curidx - 1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx - 1 - ncols, curidx - 1 - ncols - nrows*ncols);
     
     return;
} // end fcn to add four EVEN tets

void addOddTets(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     assert(i > 0 && j > 0 && k > 0);
     int curidx = i + ncols*j + ncols*nrows*k;
      
     VertexVector vertices(4);
     vertices[0] = curidx;  vertices[3] = curidx;
     vertices[1] = -1; vertices[2] = -1; 
     VertexVector::const_iterator bg = vertices.begin();
     VertexVector::const_iterator end = vertices.end();
    
     int v1, v2, v3, v4;
     double value, value2;  // max of value and value 2 is the fcn value. 

     // top vertex (i, j, k)
     v1 = curidx -1;   vertices[0] = v1;
     v2 = curidx - ncols;  vertices[1] = v2;  
     value = std::max(fcnvalues.at(v1),fcnvalues.at(v2));
     
     v3 = curidx - nrows*ncols;  vertices[2] = v3;
     v4 = curidx; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));
     
     // top vertex (i-1, j-1, k)
     v3 = curidx - 1 - ncols - nrows*ncols;  vertices[2] = v3;
     v4 = curidx - 1 -ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));

     // top vertex (i, j-1, k-1)
     v1 = curidx - nrows*ncols;   vertices[0] = v1;
     v2 = curidx - 1 - ncols - nrows*ncols;  vertices[1] = v2;  
     value = std::max(fcnvalues.at(v1),fcnvalues.at(v2));
     
     v3 = curidx - ncols;  vertices[2] = v3;
     v4 = curidx -ncols - nrows*ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
     
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));

     // top vertex (i-1, j, k-1)
     v3 = curidx - 1;  vertices[2] = v3;
     v4 = curidx -1 - nrows * ncols; vertices[3] = v4;
     value2 = std::max(fcnvalues.at(v3),fcnvalues.at(v4));
    
     filtr.push_back(Smplx(bg, bg + 4, std::max(value,value2)));
      

     return;
} // end fcn addEvenTets

void addTriNTet(Fltr& filtr, const std::vector<double> fcnvalues, 
                  const int ncols, const int nrows, int i, int j, int k)
{
     int curidx = i + ncols*j + ncols*nrows*k;
     
     // ... consider two cases for the cubical decomposition:
     if ((i+j+k)%2 == 0)
     {
	// ... EVEN BOX
        if (i > 0 && j > 0) // top
        {
     
           addTri(filtr, fcnvalues, curidx, curidx - ncols - 1, curidx - ncols);
           addTri(filtr, fcnvalues, curidx, curidx -1, curidx - ncols -1); 
        }
        if (i > 0 && k > 0) // back
        {
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols -1, curidx - 1);
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols, curidx - nrows*ncols -1);
        }
        
        if (j > 0 && k > 0) // right
        {
           addTri(filtr, fcnvalues, curidx, curidx - nrows*ncols - ncols, curidx - nrows*ncols);
           addTri(filtr, fcnvalues, curidx, curidx - ncols, curidx - nrows*ncols - ncols);
           
           if (i > 0) // middle
           {
              addTri(filtr, fcnvalues, curidx, curidx - ncols -1, curidx - ncols - nrows*ncols);
              addTri(filtr, fcnvalues, curidx, curidx -1 -nrows*ncols, curidx - ncols -1);
              addTri(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - ncols - nrows*ncols, curidx);
              addTri(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols);
             
              // ... add center tets 
              addTet(filtr, fcnvalues, curidx -1 - nrows*ncols, curidx - 1 - ncols, curidx - ncols - nrows*ncols, curidx);
              // ... add remaining tets 
              addEvenTets(filtr, fcnvalues, ncols, nrows, i, j, k);
           }
        }
     } // end if for even case 
     else {
        // ... ODD CASE
        if (i > 0 && j > 0) // top
        {
           addTri(filtr, fcnvalues, curidx -1, curidx - ncols, curidx);
           addTri(filtr, fcnvalues, curidx -1, curidx - ncols -1, curidx - ncols);
        }

        if (i > 0 && k > 0) // back
        {
           addTri(filtr, fcnvalues, curidx -1, curidx - nrows*ncols, curidx - nrows*ncols - 1);
           addTri(filtr, fcnvalues, curidx -1, curidx, curidx - nrows*ncols);
        }  
        
        if (j > 0 && k > 0) // right
        {
           addTri(filtr, fcnvalues, curidx - ncols, curidx - nrows*ncols, curidx - ncols - nrows*ncols);
           addTri(filtr, fcnvalues, curidx - ncols, curidx, curidx - nrows*ncols);
           
           if ( i > 0) // middle
           { 
              addTri(filtr, fcnvalues, curidx -1, curidx - ncols, curidx - nrows*ncols);
              addTri(filtr, fcnvalues, curidx -1, curidx - nrows*ncols - ncols - 1, curidx - ncols);
              addTri(filtr, fcnvalues, curidx - nrows*ncols, curidx - nrows*ncols - ncols -1, curidx - ncols);
              addTri(filtr, fcnvalues, curidx - nrows*ncols, curidx - 1, curidx - nrows*ncols - ncols -1);
               
              // ... add central tet
              addTet(filtr, fcnvalues, curidx -1, curidx - ncols, curidx - nrows*ncols, curidx - nrows*ncols -ncols -1);
              // ... add remaining tets
              addOddTets(filtr, fcnvalues, ncols, nrows, i, j, k);
           }
        } // end for through j k positive
     } // end else through odd case.

    return;
} // end function addTriangles


int simplicesFromGrid(Fltr& filtr, const std::string& infile)
{
  std::ifstream in(infile.c_str());
  std::string	line;
  int nrows, ncols, ndimz;
  if (std::getline(in,line))
  {
    std::stringstream linestream(line);
    if (linestream >> nrows)
    {  if (linestream >> ncols)
       { if (linestream >> ndimz)
         {; //std::cout << nrows << " rows and " << ncols << " columns." << std::endl;
         }
          else
            ndimz = 1;  // to make backwards compatible with 2d grid
       }
    }
    else
      return 1;
  }
  
  int i = 0; // indexing the columns
  int j = 0; // indexing the rows
  int k = 0; // indexing the z dimension
  std::vector<double> fcnvalues;
  //double fcnvalues [i*j]; 

  while(std::getline(in, line))
  {
    //std::vector<Vertex> currow[ncols];
    if(line[0] == '#') continue;	// comment line
    std::stringstream	linestream(line);
    double x;
    while (linestream >> x) // each line corresponds to changing i
    {
      int curidx = i + ncols*j + ncols*nrows*k;
      fcnvalues.push_back(x); // at index i + ncols*j    
      assert(fcnvalues.at(curidx) == x);
		
      // .. add the vertex 
      std::vector<Vertex> vcont;
      vcont.push_back((Vertex)(curidx));
      filtr.push_back(Smplx(vcont, fcnvalues.at(curidx))); 

      // .. NEXT, Add the edges:
      addAllEdges(filtr, fcnvalues, ncols, nrows, i, j, k);
      // ... now add the triangles:
      addTriNTet(filtr, fcnvalues, ncols, nrows, i, j, k);

      ++i; // advance column
    } // end inner while loop, which iterates through i (a row / line)

    // ... advance row / z value
    i = 0;
    ++j;
    if (j > nrows -1)
    {
	j = 0;
        ++k;
    }
  } // end while 
  in.close();
  
  return 0;
} // end simplicesFromGrid function



extern "C" {

	void grid(int* input)
	{
	#ifdef LOGGING
		//rlog::RLogInit(argc, argv);

		stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
		//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
		//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
	#endif

		// ... set up the input
		
		bool printstatus=input[0];
		
		std::string infilename;
		infilename = "inputDionysus.txt";

		Fltr f;
		simplicesFromGrid(f, infilename); // fill the simplices

		f.sort(Smplx::DataComparison()); // initialize filtration
		Persistence p(f); // initialize persistence
		p.pair_simplices(printstatus); // pair simplices
		// TODO: why doesn't this work? rLog(rlmain, "testing");   
 
		Persistence::SimplexMap<Fltr>   m = p.make_simplex_map(f);
		std::map<Dimension, PDgm> dgms;
		init_diagrams(dgms, p.begin(), p.end(), 
					  evaluate_through_map(m, Smplx::DataEvaluator()),
					  evaluate_through_map(m, Smplx::DimensionExtractor()));

		std::ofstream outfile;
		outfile.open("outputDionysus.txt");
		outfile << 0 << std::endl << dgms[0] << std::endl; // print 0-dim diagram
		outfile << 1 << std::endl << dgms[1] << std::endl; // print 1-dim diagram
		outfile << 2 << std::endl << dgms[2] << std::endl; // print 1-dim diagram
				
		// TODO: remove this line, this is just for testing
		//outfile << f;  // add the filter
	}



	void rips(int* dimInput, double* maxInput, int* printInput)
	{
		bool printstatus=printInput[0];
		Dimension               skeleton;
		DistanceType            max_distance;
		std::string             infilename, diagram_name;
		
		infilename = "inputDionysus.txt";
		diagram_name="outputDionysus.txt";
		skeleton= dimInput[0];
		max_distance= maxInput[0];
		
		std::ofstream           diagram_out(diagram_name.c_str());
		
		/*if (printstatus){
			Rprintf("Diagram %s \n", diagram_name.c_str());
			//std::cout << "Diagram:         " << diagram_name << std::endl;
		}*/
		PointContainer          points;
		read_points(infilename, points);

		PairDistances           distances(points);
		Generator               rips(distances);
		Generator::Evaluator    size(distances);
		FltrR                    f;
	
		// Generate 2-skeleton of the Rips complex for epsilon = 50
		rips.generate(skeleton, max_distance, make_push_back_functor(f));
		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", f.size());
			//std::cout << "# Generated complex of size: " << f.size() << std::endl;
		}
		// Generate filtration with respect to distance and compute its persistence
		f.sort(Generator::Comparison(distances));

		Timer persistence_timer; persistence_timer.start();
		PersistenceR p(f);
		p.pair_simplices(printstatus);
		persistence_timer.stop();

	#if 1
		// Output cycles
		PersistenceR::SimplexMap<FltrR>   m = p.make_simplex_map(f);
		for (PersistenceR::iterator cur = p.begin(); cur != p.end(); ++cur)
		{
			const PersistenceR::Cycle& cycle = cur->cycle;

			if (!cur->sign())        // only negative simplices have non-empty cycles
			{
				PersistenceR::OrderIndex birth = cur->pair;      // the cycle that cur killed was born when we added birth (another simplex)

				const SmplxR& b = m[birth];
				const SmplxR& d = m[cur];
			
				// if (b.dimension() != 1) continue;
				// std::cout << "Pair: (" << size(b) << ", " << size(d) << ")" << std::endl;
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " " << size(d) << std::endl;
			} else if (cur->unpaired())    // positive could be unpaired
			{
				const SmplxR& b = m[cur];
				// if (b.dimension() != 1) continue;
			
				// std::cout << "Unpaired birth: " << size(b) << std::endl;
				// cycle = cur->chain;      // TODO
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " inf" << std::endl;
			}

			// Iterate over the cycle
			// for (PersistenceR::Cycle::const_iterator si =  cycle.begin();
			//                                                          si != cycle.end();     ++si)
			// {
			//     const SmplxR& s = m[*si];
			//     //std::cout << s.dimension() << std::endl;
			//     const SmplxR::VertexContainer& vertices = s.vertices();          // std::vector<Vertex> where Vertex = Distances::IndexType
			//     AssertMsg(vertices.size() == s.dimension() + 1, "dimension of a simplex is one less than the number of its vertices");
			//     std::cout << vertices[0] << " " << vertices[1] << std::endl;
			// }
		}
	#endif
	
		if (printstatus){	
			persistence_timer.check("# Persistence timer");
		}

	}




	void ripsArbit(int* dimInput, double* maxInput, int* printInput)
	{
		bool printstatus=printInput[0];
		Dimension               skeleton;
		DistanceTypeA            max_distance;
		std::string             infilename, diagram_name;

		infilename = "inputDionysus.txt";
		diagram_name="outputDionysus.txt";
		skeleton= dimInput[0];
		max_distance= maxInput[0];
		
		std::ofstream           diagram_out(diagram_name.c_str());
		/*if (printstatus){
			Rprintf("Diagram %s \n", diagram_name.c_str());			
			//std::cout << "Diagram:         " << diagram_name << std::endl;
		}*/
		PointContainer          points;
		read_points2(infilename, points);

		PairDistancesA           distances(points);
		GeneratorA               rips(distances);
		GeneratorA::Evaluator    size(distances);
		FltrRA                    f;
	
		// Generate 2-skeleton of the Rips complex for epsilon = 50
		rips.generate(skeleton, max_distance, make_push_back_functor(f));
		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", f.size());
			//std::cout << "# Generated complex of size: " << f.size() << std::endl;
		}
		// Generate filtration with respect to distance and compute its persistence
		f.sort(GeneratorA::Comparison(distances));

		Timer persistence_timer; persistence_timer.start();
		PersistenceR p(f);
		p.pair_simplices(printstatus);
		persistence_timer.stop();

	#if 1
		// Output cycles
		PersistenceR::SimplexMap<FltrRA>   m = p.make_simplex_map(f);
		for (PersistenceR::iterator cur = p.begin(); cur != p.end(); ++cur)
		{
			const PersistenceR::Cycle& cycle = cur->cycle;

			if (!cur->sign())        // only negative simplices have non-empty cycles
			{
				PersistenceR::OrderIndex birth = cur->pair;      // the cycle that cur killed was born when we added birth (another simplex)

				const SmplxRA& b = m[birth];
				const SmplxRA& d = m[cur];
			
				// if (b.dimension() != 1) continue;
				// std::cout << "Pair: (" << size(b) << ", " << size(d) << ")" << std::endl;
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " " << size(d) << std::endl;
			} else if (cur->unpaired())    // positive could be unpaired
			{
				const SmplxRA& b = m[cur];
				// if (b.dimension() != 1) continue;
			
				// std::cout << "Unpaired birth: " << size(b) << std::endl;
				// cycle = cur->chain;      // TODO
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " inf" << std::endl;
			}

			// Iterate over the cycle
			// for (PersistenceR::Cycle::const_iterator si =  cycle.begin();
			//                                                          si != cycle.end();     ++si)
			// {
			//     const SmplxRA& s = m[*si];
			//     //std::cout << s.dimension() << std::endl;
			//     const SmplxR::VertexContainer& vertices = s.vertices();          // std::vector<Vertex> where Vertex = Distances::IndexType
			//     AssertMsg(vertices.size() == s.dimension() + 1, "dimension of a simplex is one less than the number of its vertices");
			//     std::cout << vertices[0] << " " << vertices[1] << std::endl;
			// }
		}
	#endif
		if (printstatus){	
			persistence_timer.check("# Persistence timer");
		}
	}



	void bottleneck(double* out_name)
	{
		// ... set up the input
		std::string filename1;
		std::string filename2;
		filename1 = "inputDionysus.txt";
		filename2 = "inputDionysus2.txt";
				
// 		po::options_description hidden("Hidden options");
// 		hidden.add_options()
// 			("input-file1",  po::value<std::string>(&filename1), "The first collection of persistence diagrams")
// 			("input-file2",  po::value<std::string>(&filename2), "The second collection of persistence diagrams");
// 
// 		po::positional_options_description p;
// 		p.add("input-file1", 1);
// 		p.add("input-file2", 2);
// 
// 		po::options_description all; all.add(hidden);
// 
// 		po::variables_map vm;
// 		po::store(po::command_line_parser(argc, argv).
// 					  options(all).positional(p).run(), vm);
// 		po::notify(vm);
// 
// 		if (!vm.count("input-file1") || !vm.count("input-file2"))
// 		{
// 			std::cout << "Usage: " << argv[0] << " input-file1 input-file2" << std::endl;
// 			return 1;
// 		}

		PDgmB dgm1, dgm2;
		read_diagram(dgm1, filename1);
		read_diagram(dgm2, filename2);
// 		std::cout << "Size dgm1: " << dgm1.size() << std::endl;
// 		std::cout << "Size dgm2: " << dgm2.size() << std::endl;
// 
// 		std::cout << "Distance: " << bottleneck_distance(dgm1, dgm2) << std::endl;
	
		out_name[0]=bottleneck_distance(dgm1, dgm2);
	}


	void wasserstein(double* inputP, double* out_name)
	{
		// ... set up the input
		std::string filename1;
		std::string filename2;
		filename1 = "inputDionysus.txt";
		filename2 = "inputDionysus2.txt";
		
		unsigned p=inputP[0];		

		PDgmB dgm1, dgm2;
		read_diagram(dgm1, filename1);
		read_diagram(dgm2, filename2);
	
		out_name[0]=wasserstein_distance(dgm1, dgm2, p);
	}


  	// KDE function on a Grid
	void kde(double *XX, int *pNN, int *pDD, double *Grid, int *pMM, double *hh, double *out){
	    double *pp= new double[pDD[0]];
		double pi=3.141593;
		double den=0.0;
		
		den=pow(hh[0], pDD[0]) * pow( 2*pi  , (pDD[0]/2.0));
		
		for (int m=1; m<=pMM[0]; m++) {
	 		for (int d=1; d<=pDD[0]; d++) {			
				pp[d-1]=ReadMat(Grid, pMM, pDD, m, d);
			}		
			out[m-1]=oneKernel(pp, XX, pNN, pDD, hh);
   			out[m-1]=out[m-1]/ den;
   		}
				
		delete[] pp;
	}


   	// kernel Dist function on a Grid
	void kdeDist(double *XX, int *pNN, int *pDD, double *Grid, int *pMM, double *hh, double *out){
	    double *pp= new double[pDD[0]];
		double first=0.0;
		double second=1.0;
	    double *third= new double[pMM[0]];
		
		for (int i=1; i<=pNN[0]; i++) {
	 		for (int d=1; d<=pDD[0]; d++) {			
				pp[d-1]=ReadMat(XX, pNN, pDD, i, d);
			}		
			first=first+ oneKernel(pp, XX, pNN, pDD, hh);
		}
		first=first/pNN[0];
		
		for (int m=1; m<=pMM[0]; m++) {
	 		for (int d=1; d<=pDD[0]; d++) {			
				pp[d-1]=ReadMat(Grid, pMM, pDD, m, d);
			}		
			third[m-1]=oneKernel(pp, XX, pNN, pDD, hh);
   			out[m-1]= std::sqrt(first+second - 2* third[m-1]  );
   		}
   		
		delete[] pp;
		delete[] third;
	}



} //end extern
