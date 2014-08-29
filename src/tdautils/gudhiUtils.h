#include <R.h>
#include <R_ext/Print.h>

#include <gudhi/io.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>

typedef int        Vertex_handle;
typedef double     Filtration_value;
