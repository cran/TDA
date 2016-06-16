#include <R.h>
#include <R_ext/Print.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>
#include <gudhi/Alpha_complex.h>

#include <utilities/timer.h>

typedef int        Vertex_handle;
typedef double     Filtration_value;



template<typename SimplexTree>
void computePersistenceGUDHI(SimplexTree& simplexTree,
		const unsigned coeffFieldCharacteristic, const double minPersistence,
		const unsigned maxdimension,
		std::vector< std::vector< std::vector< double > > > &persDgm,
		bool printProgress) {

	Timer persistence_timer;
	persistence_timer.start();

	// Compute the persistence diagram of the complex
	Gudhi::persistent_cohomology::Persistent_cohomology< Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp > pcoh(simplexTree);
	pcoh.init_coefficients(coeffFieldCharacteristic); //initializes the coefficient field for homology

	pcoh.compute_persistent_cohomology(minPersistence); //compute persistent homology

	
	std::vector< double > dgmPoint(2);
	std::vector< std::vector< double > > dgm = 
		pcoh.memory_output_diagram< std::vector< std::vector< double > > >();
	persDgm.resize(maxdimension + 1);
	for (unsigned rowIdx = 0; rowIdx < dgm.size(); ++rowIdx) {
		dgmPoint[0] = dgm[rowIdx][2];
		dgmPoint[1] = dgm[rowIdx][3];
		persDgm[dgm[rowIdx][1]].push_back(dgmPoint);
	}

	// write diagram on the output file
	//	pcoh.write_output_diagram(diagram_name);

	// or write the most persistent points in diagram (max_num_bars should be an input) 
	//  pcoh.most_persistent_bars(diagram, *max_num_bars);

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}
