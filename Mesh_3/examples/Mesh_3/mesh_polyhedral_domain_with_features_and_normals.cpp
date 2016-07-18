// #define CGAL_MESH_3_DEBUG_FACET_CRITERIA 1
#define CGAL_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/fandisk.off";
  const double angle = (argc>2)?atof(argv[2]):160;
  // Create domain
  Mesh_domain domain(fname);

  // Get sharp features
  domain.detect_features(20);

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.1,
                         facet_angle = 25);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  dump_c3t3(c3t3, "before");

  typedef Mesh_criteria::Facet_criteria::Visitor Facet_criteria_visitor;
  typedef CGAL::Mesh_3::Facet_normals_criterion<
    C3t3,
    Facet_criteria_visitor
    > Criterion;

  criteria.add_facet_criterion(new Criterion(c3t3, angle));
  CGAL::refine_mesh_3(c3t3, domain, criteria);

  dump_c3t3(c3t3, "after");
}
