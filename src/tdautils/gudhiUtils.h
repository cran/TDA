#ifndef __GUDHIUTILS_H__
#define __GUDHIUTILS_H__

// #include <R.h>
// #include <R_ext/Print.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>
#include <gudhi/Alpha_complex.h>

//for CGAL
#include <tdautils/cgalUtils.h>

#include <utilities/timer.h>

typedef int        Vertex_handle;
typedef double     Filtration_value;



// FiltrationDiag in GUDHI
/** \brief Construct the persistence diagram from the filtration using library
 *         GUDHI.
 *
 * @param[out] void             Void
 * @param[in]  smplxTree        A simplex tree
 * @param[in]  maxdimension     Max dimension of the homological features to be
 *                              computed
 * @param[in]  location         Are location of birth point, death point, and
 *                              representative cycles returned?
 * @param[in]  printProgress    Is progress printed?
 * @param[in]  persDgm          Memory space for the resulting persistence
 *                              diagram
 * @param[in]  persLoc          Memory space for the resulting birth points and
 *                              death points
 * @param[in]  persCycle        Memory space for the resulting representative
 *                              cycles
 */
template< typename SimplexTree >
inline void FiltrationDiagGudhi(
    SimplexTree    & smplxTree,
    const unsigned   coeffFieldCharacteristic,
    const double     minPersistence,
    const unsigned   maxdimension,
    bool             printProgress,
    std::vector< std::vector< std::vector< double > > > & persDgm
) {

	Timer persistence_timer;
	persistence_timer.start();

	// Compute the persistence diagram of the complex
  Gudhi::persistent_cohomology::Persistent_cohomology<
      Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp >
	    pcoh(smplxTree);
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



// RipsFiltration in GUDHI
/** \brief Construct the Rips filtration constructed on the input set of points
 *         using library GUDHI.
 *
 * @param[out] SimplexTree    A simplex tree
 * @param[in]  X              Either an nxd matrix of coordinates,
 *                            or an nxn matrix of distances of points
 * @param[in]  maxdimension   Max dimension of the homological features to be
 *                            computed.
 * @param[in]  maxscale       Threshold for the Rips complex
 * @param[in]  printProgress  Is progress printed?
 * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
 *                            diagram. Diagram must point to enough memory space for
 *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
 *                            write nothing after.
 */
template< typename SimplexTree, typename RealMatrix, typename Print >
inline SimplexTree RipsFiltrationGudhi(
    const RealMatrix & X,          //points to some memory space
    const unsigned     nSample,
    const unsigned     nDim,
    const int          maxdimension,
    const double       maxscale,
    const bool         printProgress,
    const Print      & print
) {

  // Turn the input points into a range of points
  typedef std::vector< double > Point_t;
  std::vector< Point_t > point_set =
    TdaToStl< std::vector< Point_t > >(X, nSample, nDim);


  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(point_set, maxscale
    , euclidean_distance<Point_t>);

  // Construct the Rips complex in a Simplex Tree
  Gudhi::Simplex_tree<> st;
  st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
  st.expansion(maxdimension + 1); // expand the graph until dimension dim_max

  if (printProgress) {
    print("# Generated complex of size: %d \n", st.num_simplices());
  }

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  return st;
}



//---------------------------------------------------------------------------------------------------------------------
// gudhi type definition
typedef Gudhi::Simplex_tree<>::Vertex_handle Simplex_tree_vertex;
typedef std::map<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex > Alpha_shape_simplex_tree_map;
typedef std::pair<Alpha_shape_3::Vertex_handle, Simplex_tree_vertex> Alpha_shape_simplex_tree_pair;
typedef std::vector< Simplex_tree_vertex > Simplex_tree_vector_vertex;

Vertex_list fromCell(const Cell_handle& ch)
{
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++)
  {
    the_list.push_back(ch->vertex(i));
  }
  return the_list;
}
Vertex_list fromFacet(const Facet& fct)
{
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++)
  {
    if (fct.second != i)
    {
      the_list.push_back(fct.first->vertex(i));
    }
  }
  return the_list;
}
Vertex_list fromEdge(const Edge_3& edg)
{
  Vertex_list the_list;
  for (auto i = 0; i < 4; i++)
  {
    if ((edg.second == i) || (edg.third == i))
    {
      the_list.push_back(edg.first->vertex(i));
    }
  }
  return the_list;
}
Vertex_list fromVertex(const Alpha_shape_3::Vertex_handle& vh)
{
  Vertex_list the_list;
  the_list.push_back(vh);
  return the_list;
}



// AlphaShapeFiltration in GUDHI
/** \brief Construct the alpha shape filtration constructed on the input set of
 *         points using library GUDHI.
 *
 * @param[out] SimplexTree    A simplex tree
 * @param[in]  X              An nxd matrix of coordinates
 * @param[in]  printProgress  Is progress printed?
 */
template< typename SimplexTree, typename RealMatrix, typename Print >
inline SimplexTree AlphaShapeFiltrationGudhi(
    const RealMatrix & X,
    const bool         printProgress,
    Print            & print,
    RealMatrix       & coordinates
) {

  coordinates = RealMatrix(X.nrow(), X.ncol());
  const unsigned nRow = coordinates.nrow();

  // Turn the input points into a range of points
  std::list<Point_3> lp = RcppToCGALPoint3< std::list<Point_3> >(X);


  // alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
  Alpha_shape_3 as(lp.begin(), lp.end(), 0, Alpha_shape_3::GENERAL);
  //std::cout << "Alpha shape computed in GENERAL mode" << std::endl;

  // filtration with alpha values from alpha shape
  std::vector<Object> the_objects;
  std::vector<Alpha_value_type> the_alpha_values;

  Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects), std::back_inserter(the_alpha_values));

  as.filtration_with_alpha_values(disp);
  //std::cout << "filtration_with_alpha_values returns : " << the_objects.size() << " objects" << std::endl;

  Alpha_shape_3::size_type count_vertices = 0;
  Alpha_shape_3::size_type count_edges = 0;
  Alpha_shape_3::size_type count_facets = 0;
  Alpha_shape_3::size_type count_cells = 0;

  // Loop on objects vector
  Vertex_list vertex_list;
  SimplexTree simplex_tree;
  Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
  std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
  int dim_max = 0;
  Filtration_value filtration_max = 0.0;
  for (auto object_iterator : the_objects)
  {
    // Retrieve Alpha shape vertex list from object
    if (const Cell_handle* cell = CGAL::object_cast<Cell_handle>(&object_iterator))
    {
      vertex_list = fromCell(*cell);
      count_cells++;
      if (dim_max < 3) {
        dim_max = 3; // Cell is of dim 3
      }
    }
    else if (const Facet* facet = CGAL::object_cast<Facet>(&object_iterator))
    {
      vertex_list = fromFacet(*facet);
      count_facets++;
      if (dim_max < 2) {
        dim_max = 2; // Facet is of dim 2
      }
    }
    else if (const Edge_3* edge = CGAL::object_cast<Edge_3>(&object_iterator))
    {
      vertex_list = fromEdge(*edge);
      count_edges++;
      if (dim_max < 1) {
        dim_max = 1; // Edge_3 is of dim 1
      }
    }
    else if (const Alpha_shape_3::Vertex_handle* vertex = CGAL::object_cast<Alpha_shape_3::Vertex_handle>(&object_iterator))
    {
      count_vertices++;
      vertex_list = fromVertex(*vertex);
    }
    // Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
    Simplex_tree_vector_vertex the_simplex_tree;
    for (auto the_alpha_shape_vertex : vertex_list)
    {
      Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
      if (the_map_iterator == map_cgal_simplex_tree.end())
      {
        // alpha shape not found
        Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();
        //std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] not found - insert " << vertex << std::endl;
        the_simplex_tree.push_back(vertex);
        map_cgal_simplex_tree.insert(Alpha_shape_simplex_tree_pair(the_alpha_shape_vertex, vertex));

        // added by Jisu KIM, 2018-04-23
        // extract coordinates of the corresponding vertex
        coordinates[vertex] = the_alpha_shape_vertex->point()[0];
        coordinates[vertex + nRow] = the_alpha_shape_vertex->point()[1];
        coordinates[vertex + 2 * nRow] = the_alpha_shape_vertex->point()[2];
      }
      else
      {
        // alpha shape found
        Simplex_tree_vertex vertex = the_map_iterator->second;
        //std::cout << "vertex [" << the_alpha_shape_vertex->point() << "] found in " << vertex << std::endl;
        the_simplex_tree.push_back(vertex);
      }
    }
    // Construction of the simplex_tree
    Filtration_value filtr = std::sqrt(*the_alpha_value_iterator);
    //std::cout << "filtration = " << filtr << std::endl;
    if (filtr > filtration_max) {
      filtration_max = filtr;
    }
    simplex_tree.insert_simplex(the_simplex_tree, filtr);
    if (the_alpha_value_iterator != the_alpha_values.end()) {
      ++the_alpha_value_iterator;
    }
    else {
      //std::cout << "This shall not happen" << std::endl;
    }
  }
  simplex_tree.set_filtration(filtration_max);
  simplex_tree.set_dimension(dim_max);

  if (printProgress) {
    print("# Generated complex of size: %d \n", simplex_tree.num_simplices());
  }
  //std::cout << "vertices \t\t" << count_vertices << std::endl;
  //std::cout << "edges \t\t"    << count_edges << std::endl;
  //std::cout << "facets \t\t"   << count_facets << std::endl;
  //std::cout << "cells \t\t"    << count_cells << std::endl;


  //std::cout << "Information of the Simplex Tree: " << std::endl;
  //std::cout << "  Number of vertices = " << simplex_tree.num_vertices() << " ";
  //std::cout << "  Number of simplices = " << simplex_tree.num_simplices() << std::endl << std::endl;
  //std::cout << "  Dimension = " << simplex_tree.dimension() << " ";
  //std::cout << "  filtration = " << simplex_tree.filtration() << std::endl << std::endl;
  //std::cout << "Iterator on vertices: " << std::endl;
  //for( auto vertex : simplex_tree.complex_vertex_range() )
  //{ std::cout << vertex << " "; }

  // Sort the simplices in the order of the filtration
  simplex_tree.initialize_filtration();

  //std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;

  return simplex_tree;
}



// AlphaComplexFiltration in GUDHI
/** \brief Construct the alpha complex filtration constructed on the input set
 *         of points using library GUDHI.
 *
 * @param[out] SimplexTree    A simplex tree
 * @param[in]  X              An nxd matrix of coordinates
 * @param[in]  printProgress  Is progress printed?
 */
template< typename SimplexTree, typename RealMatrix, typename Print >
inline SimplexTree AlphaComplexFiltrationGudhi(
    const RealMatrix & X,
    const bool         printProgress,
    Print            & print
) {

	using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>;
	using Point = Kernel::Point_d;

	// Turn the input points into a range of points
	std::list<Point> lp = RcppToCGALPointD< std::list< Point > >(X);
	

	//Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_points(lp, maxalphasquare);
	Gudhi::alphacomplex::Alpha_complex<Kernel>
		alpha_complex_from_points(lp, std::numeric_limits<double>::infinity());

	if (printProgress) {
		print("# Generated complex of size: %d \n", alpha_complex_from_points.num_simplices());
	}

	// Sort the simplices in the order of the filtration
	alpha_complex_from_points.initialize_filtration();

	return alpha_complex_from_points;
}



# endif // __GUDHIUTILS_H__
