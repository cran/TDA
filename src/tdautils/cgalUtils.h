#include <R.h>
#include <R_ext/Print.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/iterator.h>

// Alpha_shape_3 templates type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Alpha_shape_vertex_base_3<Kernel>             Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel>               Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>         Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds>          Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>                Alpha_shape_3;

// From file type definition
typedef Kernel::Point_3                                     Point_3;

// filtration with alpha values needed type definition
typedef Alpha_shape_3::FT Alpha_value_type;
typedef CGAL::Object      Object;
typedef CGAL::Dispatch_output_iterator<
  CGAL::cpp11::tuple<Object, Alpha_value_type>,
  CGAL::cpp11::tuple<std::back_insert_iterator< std::vector<Object> >, std::back_insert_iterator< std::vector<Alpha_value_type> >
                     > > Dispatch;
typedef Alpha_shape_3::Cell_handle   Cell_handle;
typedef Alpha_shape_3::Facet         Facet;
typedef Alpha_shape_3::Edge          Edge_3;
typedef std::list<Alpha_shape_3::Vertex_handle> Vertex_list;
