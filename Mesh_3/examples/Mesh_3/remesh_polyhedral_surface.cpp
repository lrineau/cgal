#define CGAL_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>

// Domain 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#if CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char** argv)
{
  // Load OBJ
  std::vector<K::Point_3> points;
  std::vector<std::vector<std::size_t> > faces;

  std::ifstream in(argv[1]);
  if(!in || !CGAL::read_OBJ(in,points,faces))
    return 1;

  Polyhedron poly;
  namespace PMP = CGAL::Polygon_mesh_processing;
  PMP::orient_polygon_soup(points,faces);
  PMP::polygon_soup_to_polygon_mesh(points, faces, poly);
  if (!CGAL::is_triangle_mesh(poly)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create a vector with only one element: the pointer to the polyhedron.
  std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);

  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
  
  // Mesh criteria
  Mesh_criteria criteria(facet_angle = 25,
                         facet_distance = 0.05);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output the facets of the c3t3 to an OFF file. The facets will not be
  // oriented.
  Polyhedron out_surface_mesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_surface_mesh);

  std::ofstream off_file("out.off");
  off_file.precision(17);
  off_file << out_surface_mesh;

  return off_file.fail() ? EXIT_FAILURE : EXIT_SUCCESS;
}
