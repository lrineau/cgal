#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/SMDS_3/Dump_c3t3.h>

#include "read_polylines.h"

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::FT FT;

typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;

// Implicit Domain
typedef CGAL::Labeled_mesh_domain_3<Kernel> Implicit_domain;

// Polyhedral Domain
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Polyhedron_domain;

using Subdomain_index = int;
using Surface_patch_index = unsigned char;
using Curve_index = char;
using Corner_index = short;
using Cb = CGAL::Simplicial_mesh_cell_base_3<Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<Kernel, Subdomain_index, Surface_patch_index,
                                               Curve_index, Corner_index>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrency_tag>;
using Triangulation = CGAL::Triangulation_3<Kernel, Tds>;

class Hybrid_domain {
	const Implicit_domain& implicit_domain;
	const Polyhedron_domain& polyhedron_domain;

	public:
	Hybrid_domain(const Implicit_domain& implicit_domain, const Polyhedron_domain& polyhedron_domain)
	              : implicit_domain(implicit_domain) , polyhedron_domain(polyhedron_domain) {}
	
	// types required by the 'MeshDomain_3' concept
	typedef int Surface_patch_index;
	typedef int Subdomain_index;
	typedef int Index;
	typedef Kernel R;
	typedef Kernel::Point_3 Point_3;
	typedef std::tuple<Point_3, Index, int> Intersection;
  
	CGAL::Bbox_3 bbox() const { return implicit_domain.bbox() + polyhedron_domain.bbox(); }

	struct Construct_initial_points {
		Construct_initial_points(const Hybrid_domain& domain) : r_domain_(domain) {}
		template<class OutputIterator> OutputIterator operator()(OutputIterator pts, const int n = 20) const {
			//construct initial points on implicit domain
			typedef Implicit_domain::Index Implicit_Index;
			std::vector<std::pair<Point_3, Implicit_Index> > implicit_points_vector;
			Implicit_domain::Construct_initial_points cstr_implicit_initial_points =
				r_domain_.implicit_domain.construct_initial_points_object();
			cstr_implicit_initial_points(std::back_inserter(implicit_points_vector), n/2);
			for(std::size_t i = 0, end = implicit_points_vector.size(); i<end; ++i) {
				*pts++ = std::make_pair(implicit_points_vector[i].first, 2);
			}
			//construct initial points on polyhedral domain
			typedef Polyhedron_domain::Index Polyhedron_Index;
			std::vector<std::pair<Point_3, Polyhedron_Index> > polyhedron_points_vector;
			Polyhedron_domain::Construct_initial_points cstr_polyhedron_initial_points =
				r_domain_.polyhedron_domain.construct_initial_points_object();
			cstr_polyhedron_initial_points(std::back_inserter(polyhedron_points_vector), n/2);
			for(std::size_t i = 0, end = polyhedron_points_vector.size(); i<end; ++i) {
				*pts++ = std::make_pair(polyhedron_points_vector[i].first, 1);
			}
			return pts;
		}

	private:
	const Hybrid_domain& r_domain_;
	}; // end Construct_initial_points_object

	Construct_initial_points construct_initial_points_object() const {
		return Construct_initial_points(*this);
	}

	struct Is_in_domain {
		Is_in_domain(const Hybrid_domain& domain) : r_domain_(domain) {}
		boost::optional<Subdomain_index> operator()(const Kernel::Point_3& p) const {
			boost::optional<Subdomain_index> implicit_subdomain_index = 
				r_domain_.implicit_domain.is_in_domain_object()(p);
			boost::optional<Subdomain_index> polyhedron_subdomain_index = 
				r_domain_.polyhedron_domain.is_in_domain_object()(p);

			if(!implicit_subdomain_index && polyhedron_subdomain_index)
				return 2;
			else
				return polyhedron_subdomain_index;

		}

		private:
		const Hybrid_domain& r_domain_;
	}; // end Is_in_domain_object

	Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }
  
	struct Construct_intersection {
		Construct_intersection(const Hybrid_domain& domain) : r_domain_(domain) {}
		template <typename Query> Intersection operator()(const Query& query) const {
			
			using boost::get;
		
			// intersection with implicit domain (within the polyhedral domain)
			Implicit_domain::Intersection implicit_inter =
				r_domain_.implicit_domain.construct_intersection_object()(query);
			// if found, return it
			if ( get<2>(implicit_inter) != 0 ) {
				const Point_3 inter_point = get<0>(implicit_inter);
				if ( r_domain_.polyhedron_domain.is_in_domain_object()(inter_point) ) // Inner implicit surface
					return Intersection(inter_point, 4, get<2>(implicit_inter));
			}
			
			// intersection with polyhedral domain
			Polyhedron_domain::Intersection polyhedron_inter =
				r_domain_.polyhedron_domain.construct_intersection_object()(query);
			// if found, return it
			if ( get<2>(polyhedron_inter) != 0 ) {
				const Point_3 inter_point = get<0>(polyhedron_inter);
				if ( r_domain_.implicit_domain.is_in_domain_object()(inter_point) ) // Scaffold
					return Intersection(inter_point, 4, get<2>(polyhedron_inter));
				else // Void
					return Intersection(inter_point, 5, get<2>(polyhedron_inter));
			}

			//no intersection found
			return Intersection();
		}

		private:
		const Hybrid_domain& r_domain_;
	}; // end Construct_intersection_object

	Construct_intersection construct_intersection_object() const {
		return Construct_intersection(*this);
	}

	//Index types converters
	Index index_from_surface_patch_index(const Surface_patch_index& index) const
		{ return index; }
	Index index_from_subdomain_index(const Subdomain_index& index) const
		{ return index; }
	Surface_patch_index surface_patch_index(const Index& index) const
		{ return index; }
	Subdomain_index subdomain_index(const Index& index) const
		{ return index; }
}; // end Hybrid_domain class

typedef CGAL::Mesh_domain_with_polyline_features_3<Hybrid_domain> H_domain;
typedef CGAL::Mesh_triangulation_3<H_domain, CGAL::Default, Concurrency_tag>::type H_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<H_Tr> C3t3;
typedef CGAL::Mesh_criteria_3<H_Tr>     H_Mesh_criteria;
typedef H_Mesh_criteria::Edge_criteria  H_Edge_criteria;
typedef H_Mesh_criteria::Facet_criteria H_Facet_criteria;
typedef H_Mesh_criteria::Cell_criteria  H_Cell_criteria;

FT my_implicit_surface (const Point_3& p) {
	const double x2 = p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
	double a = (x2+y2+z2-0.8); // Outer shape
	
	double scaling = (2*3.14159265359)/1.5;
	double b = sin(scaling*p.x())*cos(scaling*p.y()) +
	    sin(scaling*p.y())*cos(scaling*p.z()) +
	    sin(scaling*p.z())*cos(scaling*p.x()) -
		 0;
				
	if (a > b) {return a;} else {return b;}
}

int main() {
  const std::string fname = CGAL::data_file_path("cube.off");
  
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.bad()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  // Domains
  Polyhedron_domain polyhedron_domain(polyhedron);
  Implicit_domain sphere_domain =
    Implicit_domain::create_implicit_mesh_domain(my_implicit_surface,
                                                 Kernel::Sphere_3(Point_3(0, 0, 0), FT(9)));
  H_domain domain(sphere_domain, polyhedron_domain);
  
  // Polyline features
  const char* lines_fname = "polylines.txt";
  std::vector<std::vector<Point_3> > featured_curves;
  if (!read_polylines(lines_fname, featured_curves)) {
    std::cerr << "Error: Cannot read file " << lines_fname << std::endl;
    return EXIT_FAILURE;
  }
 
  // Add features for protection
  domain.add_features(featured_curves.begin(), featured_curves.end());

  // Criteria
  H_Edge_criteria edge_criteria(0.05);
  H_Facet_criteria facet_criteria(30, 0.05, 0.025); // angle, size, approximation
  H_Cell_criteria cell_criteria(2, 0.05); // radius-edge ratio, size
  H_Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);
  
  // Mesh generation (without optimization)
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_perturb().no_exude());
  // Output
  dump_c3t3(c3t3, "out");
  
  return EXIT_SUCCESS;
}