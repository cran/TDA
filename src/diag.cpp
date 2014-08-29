#include <R.h>
#include <R_ext/Print.h>


//for Rips
#include <tdautils/ripsL2.h>
#include <tdautils/ripsArbit.h>

//for grid
#include <tdautils/gridUtils.h>

//for kernel density
#include <tdautils/kernelUtils.h>

//for bottleneck and Wasserstein
#include <tdautils/bottleneckUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

// for phat
#include <tdautils/phatUtils.h>



extern "C" {



	// grid function by Brittany T. Fasy
	// modified by Jisu Kim for arbitrary dimension & using memory as an input & setting maximum dimension
	void grid(double *FUNvaluesInput, int *gridDimensionInput, int *gridNumberInput, int *maxdimensionInput, char** decompositionInput, char** libraryInput, int* printInput)
	{
	#ifdef LOGGING
		//rlog::RLogInit(argc, argv);

		stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
		//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
		//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
	#endif

		bool printstatus=printInput[0];
		Fltr f;

		// Generate simplicial complex from function values and grid
		if (decompositionInput[0][0] == '5')
		{
			simplicesFromGrid(f, FUNvaluesInput, gridNumber( gridDimensionInput, gridNumberInput ), (*maxdimensionInput)+1 ); // fill the simplices
		}
		if (decompositionInput[0][0] == 'b')
		{
			simplicesFromGridBarycenter(f, FUNvaluesInput, gridNumber( gridDimensionInput, gridNumberInput ), (*maxdimensionInput)+1 ); // fill the simplices
		}
		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", f.size());
		}


		// Sort simplicial complex
		f.sort(Smplx::DataComparison()); // initialize filtration


		// Compute persistent homology from sorted simplicial complex
		std::map<Dimension, PDgm> dgms;
		std::vector< std::vector< std::vector< double > > > persDgm;
		if (libraryInput[0][0] == 'D')
		{

			Persistence p(f); // initialize persistence
			p.pair_simplices(printstatus); // pair simplices
			// TODO: why doesn't this work? rLog(rlmain, "testing");   
	 
			Persistence::SimplexMap<Fltr> m = p.make_simplex_map(f);
			init_diagrams(dgms, p.begin(), p.end(), 
						  evaluate_through_map(m, Smplx::DataEvaluator()),
						  evaluate_through_map(m, Smplx::DimensionExtractor()));

			//std::vector< double > persDgmPoint(3);
			//std::map<Dimension, PDgm>::const_iterator dgmItr;
			//for (dgmItr = dgms.begin(); dgmItr != dgms.end(); ++dgmItr)
			//{
			//		if (dgmItr->first > *maxdimensionInput)
			//				break;
			//		PDgm::const_iterator dgmPtItr;
			//		for (dgmPtItr = (dgmItr->second).begin(); dgmPtItr != (dgmItr->second).end(); ++dgmPtItr)
			//		{
			//				persDgmPoint[0] = dgmItr->first;
			//				persDgmPoint[1] = dgmPtItr->x();
			//				persDgmPoint[2] = dgmPtItr->y();
			//				persDgm.push_back( persDgmPoint );
			//		}
			//}
		}
		if (libraryInput[0][0] == 'P')
		{
			persDgm = computePersistentPairsPhat(f, *maxdimensionInput);
		}


		// Output persistent diagram
		std::ofstream outfile;
		outfile.open("outputTDA.txt");
		for (int dgmsIdx=0; dgmsIdx<= *maxdimensionInput; ++dgmsIdx)
		{
			outfile << dgmsIdx << std::endl;
			if (libraryInput[0][0] == 'D')
			{
				outfile << dgms[dgmsIdx] << std::endl; // print i-dim diagram
			}
			if (libraryInput[0][0] == 'P')
			{
				std::vector< std::vector< double > >::const_iterator persdgmIdx;
				for (persdgmIdx = persDgm[ dgmsIdx ].begin(); persdgmIdx != persDgm[ dgmsIdx ].end(); ++persdgmIdx)
				{
					outfile << (*persdgmIdx)[0] << " " << (*persdgmIdx)[1] << std::endl;
				}
			}
		}
		outfile.close();
	}



	void rips(int* dimInput, double* maxInput, int* printInput)
	{
		bool printstatus=printInput[0];
		Dimension               skeleton;
		DistanceType            max_distance;
		std::string             infilename, diagram_name;
		
		infilename = "inputTDA.txt";
		diagram_name="outputTDA.txt";
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

		infilename = "inputTDA.txt";
		diagram_name="outputTDA.txt";
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



	void bottleneck(double* points1Input, int* points1NumberInput, double* points2Input, int* points2NumberInput, double* out_name)
	{
		// ... set up the input
		PDgmB dgm1, dgm2;
		read_diagram(dgm1, points1Input, points1NumberInput);
		read_diagram(dgm2, points2Input, points2NumberInput);

		out_name[0]=bottleneck_distance(dgm1, dgm2);
	}


	void wasserstein(double* points1Input, int* points1NumberInput, double* points2Input, int* points2NumberInput, int* inputP, double* out_name)
	{
		// ... set up the input
		int p=inputP[0];		

		PDgmB dgm1, dgm2;
		read_diagram(dgm1, points1Input, points1NumberInput);
		read_diagram(dgm2, points2Input, points2NumberInput);
	
		out_name[0]=wasserstein_distance(dgm1, dgm2, p);
	}


  	// KDE function on a Grid
	void kde(double *XX, int *pNN, int *pDD, double *Grid, int *pMM, double *hh, int *printProgress, double *out){
	    double *pp= new double[pDD[0]];
		double pi=3.141593;
		double den=0.0;

		int counter=0;
		int totalCount=pMM[0];
		int percentageFloor=0;	
		int tmp;

		den=pow(hh[0], pDD[0]) * pow( 2*pi  , (pDD[0]/2.0));
		
		if (printProgress[0])
		{
			Rprintf("0   10   20   30   40   50   60   70   80   90   100");
			Rprintf("\n");
			Rprintf("|----|----|----|----|----|----|----|----|----|----|\n");
			Rprintf("*");		


			for (int m=1; m<=pMM[0]; m++) {
				for (int d=1; d<=pDD[0]; d++) {			
					pp[d-1]=ReadMat(Grid, pMM, pDD, m, d);
				}		
				out[m-1]=oneKernel(pp, XX, pNN, pDD, hh);
				out[m-1]=out[m-1]/ den;
						
				//printProgress
				counter++;
				tmp=std::floor((100*counter/totalCount-percentageFloor)/2);
				if (tmp>0)
				{
					for (int aa=1; aa<=tmp; aa++)
					{
						Rprintf("*");
						percentageFloor=percentageFloor+2;
					}							
				}					   		
			}
		} else //no printProgress
		{
			for (int m=1; m<=pMM[0]; m++) {
				for (int d=1; d<=pDD[0]; d++) {			
					pp[d-1]=ReadMat(Grid, pMM, pDD, m, d);
				}		
				out[m-1]=oneKernel(pp, XX, pNN, pDD, hh);
				out[m-1]=out[m-1]/ den;					   		
			}

		}
   		
   		if (printProgress[0]) Rprintf("\n");					
		delete[] pp;
	}


   	// kernel Dist function on a Grid
	void kdeDist(double *XX, int *pNN, int *pDD, double *Grid, int *pMM, double *hh, int *printProgress, double *out){
	    double *pp= new double[pDD[0]];
		double first=0.0;
		double second=1.0;
	    double *third= new double[pMM[0]];
		
		int counter=0;
		int totalCount=pNN[0]+pMM[0];
		int percentageFloor=0;	
		int tmp;
		
		if (printProgress[0])
		{
			Rprintf("0   10   20   30   40   50   60   70   80   90   100");
			Rprintf("\n");
			Rprintf("|----|----|----|----|----|----|----|----|----|----|\n");
			Rprintf("*");		
					
			for (int i=1; i<=pNN[0]; i++) {
				for (int d=1; d<=pDD[0]; d++) {			
					pp[d-1]=ReadMat(XX, pNN, pDD, i, d);
				}		
				first=first+ oneKernel(pp, XX, pNN, pDD, hh);
		
				// printProgress
				counter++;
				tmp=std::floor((100*counter/totalCount-percentageFloor)/2);
				
				if (tmp>0)
				{
					for (int aa=1; aa<=tmp; aa++)
					{
						Rprintf("*");
						percentageFloor=percentageFloor+2;
					}							
				}				
			}
			first=first/pNN[0];
		
			for (int m=1; m<=pMM[0]; m++) {
				for (int d=1; d<=pDD[0]; d++) {			
					pp[d-1]=ReadMat(Grid, pMM, pDD, m, d);
				}		
				third[m-1]=oneKernel(pp, XX, pNN, pDD, hh);
				out[m-1]= std::sqrt(first+second - 2* third[m-1]  );

				// printProgress
				counter++;
				tmp=std::floor((100*counter/totalCount-percentageFloor)/2);
				if (tmp>0)
				{
					for (int aa=1; aa<=tmp; aa++)
					{
						Rprintf("*");
						percentageFloor=percentageFloor+2;
					}							
				}
			}   		
			Rprintf("\n");		
   		} else //no printProgress
		{
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
   		}
		delete[] pp;
		delete[] third;
	}



	// GUDHI RIPS
	/** \brief Interface for R code, construct the persistence diagram 
	  * of the Rips complex constructed on the input set of points.
	  *
	  * @param[out] void            every function called by R must return void
	  * @param[in]  point           pointer toward the coordinates of all points. Format 
	  *                             must be X11 X12 ... X1d X21 X22 ... X2d X31 ...
	  * @param[in]  dim             embedding dimension 
	  * @param[in]  num_points      number of points. The input point * must be a 
	  *                             pointer toward num_points*dim double exactly.
	  * @param[in]  rips_threshold  threshold for the Rips complex
	  * @param[in]  max_complex_dim maximal dimension of the Rips complex
	  * @param[in]  diagram         where to output the diagram. The format must be dimension birth death.
	  * @param[in]  max_num_bars    write the max_num_pairs most persistent pairs of the
	  *                             diagram. diagram must point to enough memory space for
	  *                             3*max_num_pairs double. If there is not enough pairs in the diagram,
	  *                             write nothing after.
	  */
	void 
	rips_persistence_diagram_GUDHI( double * points          //points to some memory space
							, int    * dim
							, int    * num_points
							, double * rips_threshold
							, int    * max_complex_dim
							, double * diagram         //points to some memory space
							, int    * printInput)

	{
	  bool printstatus=printInput[0];  
	  std::string diagram_name="outputTDA.txt";


	  // Turn the input points into a range of points
	  typedef std::vector<double> Point_t;
	  std::vector< Point_t > point_set(*num_points, Point_t(*dim));
	  double * curr_coordinate = points;
	  for(int i = 0; i < *num_points; ++i) {
		for(int j = 0; j < *dim; ++j) {
		  point_set[i][j] = *curr_coordinate;
		  ++curr_coordinate;
		}
	  }
	// Compute the proximity graph of the points
	  Graph_t prox_graph = compute_proximity_graph( point_set, *rips_threshold
												  , euclidean_distance<Point_t> );

	// Construct the Rips complex in a Simplex Tree
	// Construct the Rips complex in a Simplex Tree
	  Simplex_tree<> st;        
	  st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
	  st.expansion( *max_complex_dim ); // expand the graph until dimension dim_max

		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", st.num_simplices());
			// std::cout << st.num_simplices() << " simplices \n";
		}

	// Sort the simplices in the order of the filtration
	  st.initialize_filtration();

	// Compute the persistence diagram of the complex
	  int p = 2; //characteristic of the coefficient field for homology
	  double min_persistence = 0; //minimal length for persistent intervals
	  Persistent_cohomology< Simplex_tree<>, Field_Zp > pcoh( st );
	  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
	  pcoh.compute_persistent_cohomology( min_persistence ); //compute persistent homology

	// write diagram on the output file
	  pcoh.write_output_diagram(diagram_name);
	
	// or write the most persistent points in diagram (max_num_bars should be an input) 
	//  pcoh.most_persistent_bars(diagram, *max_num_bars);
	}



} //end extern
