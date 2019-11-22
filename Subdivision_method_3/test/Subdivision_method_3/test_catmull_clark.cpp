/*
 * test_catmull_clark.cpp
 *
 *  Created on: 8 Nov 2019
 *      Author: xxiao
 */

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/circulator.h>

#include <boost/lexical_cast.hpp>

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/if.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <eigen3/Eigen/Dense>

#include <dtkAbstractRationalBezierSurfaceData.h>
#include <dtkRationalBezierSurface.h>
#include <dtkContinuousGeometry.h>

#include <iterator>
#include <list>
#include <vector>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>    PolygonMesh;

typedef typename boost::graph_traits<PolygonMesh> graph_traits;

typedef typename graph_traits::vertex_descriptor vertex_descriptor;
typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
typedef typename graph_traits::edge_descriptor edge_descriptor;
typedef typename graph_traits::face_descriptor face_descriptor;

using namespace std;
using namespace CGAL;

// declaration of auxiliary functions
void onering_vertices_around_vertex(std::vector< vertex_descriptor >& onering_vertices_e,
                                    std::vector< vertex_descriptor >& onering_vertices_f,
                                    const vertex_descriptor& vtx,
                                    const PolygonMesh& pmesh);

void onering_vertices_around_face(std::vector< vertex_descriptor >& onering_vertices,
                                  const face_descriptor& f,
                                  const PolygonMesh& pmesh);

void collect_subdivision_map(std::map< face_descriptor, std::vector< halfedge_descriptor > >& subdiv_map,
                             const PolygonMesh& pmesh,
                             const PolygonMesh& pmesh_subdiv);

void index_reorder_map(std::vector< std::pair<int, int> >& indices, const int& id1, const int& id2);

void vertices_reorder_map(boost::numeric::ublas::bounded_matrix< vertex_descriptor, 5, 5 >& vtx_map,
                          const face_descriptor& f, const PolygonMesh& pmesh1, const PolygonMesh& pmesh2,
                          std::map< face_descriptor, std::vector< halfedge_descriptor > > subdiv_map_1,
                          std::map< face_descriptor, std::vector< halfedge_descriptor > > subdiv_map_2);

void vertices_limit_map(boost::numeric::ublas::bounded_matrix< Kernel::Point_3, 5, 5 >& vtx_limit_map,
                        boost::numeric::ublas::bounded_matrix< vertex_descriptor, 5, 5 > vtx_map,
                        const PolygonMesh& pmesh2);

bool check_regular_patch(const face_descriptor& f, const PolygonMesh& pmesh);

void collect_edge_vertices_LS(std::vector< std::pair<int, int> >& indices, const int& edge_id);

void edge_least_square_fitting(std::vector< Kernel::Point_3 >& bezier_points, const std::vector< Kernel::Point_3 >& edge_limit_points);

void collect_face_vertices_LS(std::vector< std::pair<int, int> >& indices);

void face_least_square_fitting(std::map<int, Kernel::Point_3 >& bezier_points, const std::vector< Kernel::Point_3 >& face_limit_points);

void bspline_to_bezier(std::map<int, Kernel::Point_3>& bezier_points, const Eigen::MatrixXd& cps_spline);

void create_bezier_control_points(std::map< std::pair<int, int>, Kernel::Point_3 >& cps_bezier, std::map<int, Kernel::Point_3> bezier_patch);

Kernel::Point_3 evaluate_bezier_surface(const double& u, const double& v,
                                        std::map< std::pair<int, int>, Kernel::Point_3 > cps_bezier);

int main(int argc, char** argv) {
  if (argc > 4) {
    cerr << "Usage: CatmullClark_subdivision [d] [filename_in] [filename_out] \n";
    cerr << "         d -- the depth of the subdivision (default: 1) \n";
    cerr << "         filename_in -- the input mesh (.off) (default: data/quint_tris.off) \n";
    cerr << "         filename_out -- the output mesh (.off) (default: result.off)" << endl;
    return 1;
  }

  int d = (argc > 1) ? boost::lexical_cast<int>(argv[1]) : 1;
  const char* in_file = (argc > 2) ? argv[2] : "data/quint_tris.off";
  const char* out_file = (argc > 3) ? argv[3] : "result.off";

  PolygonMesh pmesh;
  std::ifstream in(in_file);
  if(in.fail()) {
    std::cerr << "Could not open input file " << in_file << std::endl;
    return 1;
  }
  in >> pmesh;

  // test of surface mesh usage
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef typename graph_traits::halfedge_iterator halfedge_iterator;
  typedef typename graph_traits::face_iterator face_iterator;

  typedef CGAL::Halfedge_around_face_circulator<PolygonMesh> Halfedge_around_face_circulator;
  typedef CGAL::Vertex_around_face_circulator<PolygonMesh> Vertex_around_face_circulator;

  // mesh information
  typename graph_traits::vertices_size_type num_vertices = CGAL::num_vertices(pmesh);
  typename graph_traits::halfedges_size_type num_edges = CGAL::num_halfedges(pmesh)/2;
  typename graph_traits::faces_size_type num_facets = CGAL::num_faces(pmesh);

  std::cout << "The polygon mesh has ";
  std::cout << num_vertices << " vertices, ";
  std::cout << num_edges << " edges, ";
  std::cout << num_facets << " faces." << std::endl;

  // find subdivided vertices of a facet
  // define an auxiliary mesh to subdivide
  PolygonMesh pmesh1;
  //CGAL::copy_face_graph(pmesh, pmesh1);
  in.clear();
  in.seekg(0, ios::beg);
  in >> pmesh1;

  Subdivision_method_3::CatmullClark_subdivision(pmesh1,
      Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, pmesh1)).number_of_iterations(1));

  //===========================================================================
  // using halfedge information before and after subdivision
  //===========================================================================
  // loop over facets of the original mesh
  std::map<face_descriptor, std::vector<halfedge_descriptor> > subdiv_map_1;

  collect_subdivision_map(subdiv_map_1, pmesh, pmesh1);
/*
  for (face_iterator fitr = pmesh.faces_begin(); fitr != pmesh.faces_end(); ++fitr) {
    std::cout << "original facet " << *fitr << ": ";
    std::vector<halfedge_descriptor> halfedges = subdiv_map_1[*fitr];
    for (int j = 0; j < halfedges.size(); ++j) {
      halfedge_descriptor he = halfedges[j];
      face_descriptor f = pmesh1.face(he);
      std::cout << f << ", ";
    }
    std::cout << std::endl;
  }
*/
  // define another auxiliary mesh to subdivide
  PolygonMesh pmesh2;
  in.clear();
  in.seekg(0, ios::beg);
  in >> pmesh2;

  Subdivision_method_3::CatmullClark_subdivision(pmesh2,
        Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, pmesh2)).number_of_iterations(2));

  // loop over facets of the subdivided mesh pmesh1
  std::map<face_descriptor, std::vector<halfedge_descriptor> > subdiv_map_2;

  collect_subdivision_map(subdiv_map_2, pmesh1, pmesh2);

  // loop over original mesh and locate subdivided vertices
  for (face_iterator fitr = pmesh.faces_begin(); fitr != pmesh.faces_end(); ++fitr) {
    face_descriptor f = *fitr;

    //std::cout << "original facet " << f << ": " << std::endl;

    // rearrange the subdivided vertices in a 5 x 5 map
    boost::numeric::ublas::bounded_matrix< vertex_descriptor, 5, 5 > vtx_map; // it contains 5 x 5 vertices after two subdivisions of a facet
    vertices_reorder_map(vtx_map, f, pmesh1, pmesh2, subdiv_map_1, subdiv_map_2);

    // compute limit positions of vertices
    boost::numeric::ublas::bounded_matrix< Kernel::Point_3, 5, 5 > vtx_limit_map;
    vertices_limit_map(vtx_limit_map, vtx_map, pmesh2);

    // create Bezier patch
    bool is_regular = check_regular_patch(f, pmesh);

    typename PolygonMesh::size_type degree_face = pmesh.degree(f);

    if (is_regular == true) {
      // regular patch
      std::vector<vertex_descriptor> onering_vertices_face(16);
      onering_vertices_around_face(onering_vertices_face, f, pmesh);

      Eigen::MatrixXd cps_spline;
      cps_spline.resize(16, 3);

      for (int i = 0; i < 16; ++i) {
        vertex_descriptor vtx = onering_vertices_face[i];
        Kernel::Point_3 point = pmesh.point(vtx);
        for (int j = 0; j < 3; ++j) {
          cps_spline(i, j) = point[j];
        }
      }

      std::map<int, Kernel::Point_3> bezier_patch_regular;

      bspline_to_bezier(bezier_patch_regular, cps_spline);
/*
      // output regular bezier control points
      for (int i = 0; i < bezier_patch_regular.size(); ++i) {
        Kernel::Point_3 cp = bezier_patch_regular[i];
        std::cout << cp << std::endl;
      }
*/

      // rearrange Bezier control points
      std::map< std::pair<int, int>, Kernel::Point_3 > cps_bezier_regular;
      create_bezier_control_points(cps_bezier_regular, bezier_patch_regular);
/*
      for (int j = 0; j < 11; ++j) {
        double v = 0. + 0.1*j;
        for (int i = 0; i < 11; ++i) {
          double u = 0. + 0.1*i;
          Kernel::Point_3 point = evaluate_bezier_surface(u, v, cps_bezier_regular);
          std::cout << point << std::endl;
        }
      }
*/
      dtkContinuousGeometry::initialize();
      
      dtkAbstractRationalBezierSurfaceData* bezier_surface_data = dtkContinuousGeometry::abstractRationalBezierSurfaceData::pluginFactory().create("dtkRationalBezierSurfaceDataOn");
      // if(!bezier_surface_data) {
      //   std::cerr << "nurbs surface data could not be instantiated" << std::endl;
      //   dtkFatal() << "nurbs surface data could not be instantiated";
      // }

      dtkRationalBezierSurface nurbs_surface(bezier_surface_data);

      //bi-degree 1,2
      std::size_t dim = 3;
      std::size_t order_u = 4;
      std::size_t order_v = 4;

      double* cps = new double[(dim + 1) * order_u * order_v];

      for (int i = 0; i < bezier_patch_regular.size(); ++i) {
        Kernel::Point_3 coord = bezier_patch_regular[i];
        for (int j = 0; j < 3; ++j) {
          cps[(dim + 1)*i + j] = coord[j];
        }
        cps[(dim + 1)*i + 3] = 1.;
      }

      nurbs_surface.create(dim, order_u, order_v, cps);

      double* point = new double[3];
      for (double i = 0.; i < 1.; i+=0.01){
        for (double j = 0.; j < 1.; j+=0.01){
          nurbs_surface.evaluatePoint(i, j, point);
          Kernel::Point_3 point_1 = evaluate_bezier_surface(j, i, cps_bezier_regular);
          std::cerr << "dtk: " << point[0] << " " << point[1] << " " << point[2] << std::endl;
          //surface_points_file << point[0] << " " << point[1] << " " << point[2] << std::endl;

          std::cerr << "test: " << point_1[0] << " " << point_1[1] << " " << point_1[2] << std::endl;
          std::cout << std::endl;
        }
      }
    }
    else {
      // irregular patch, need least-square fitting approximation
      std::map<int, Kernel::Point_3> bezier_patch_ev;

      // corners use limit positions directly
      bezier_patch_ev[0] = vtx_limit_map(0, 0);
      bezier_patch_ev[3] = vtx_limit_map(0, 4);
      bezier_patch_ev[12] = vtx_limit_map(4, 0);
      bezier_patch_ev[15] = vtx_limit_map(4, 4);

      // internal edges
      // start from the edge with the smallest target id
      halfedge_descriptor he_min = pmesh.halfedge(f);
      halfedge_descriptor he_tmp = pmesh.halfedge(f);
      for (int i = 0; i < degree_face; ++i) {
        if ( CGAL::target(he_tmp, pmesh) < CGAL::target(he_min, pmesh) ) {
          he_min = he_tmp;
        }
        he_tmp = CGAL::next(he_tmp, pmesh);
      }

      halfedge_descriptor he = he_min;
      // loop over the halfedges of the facet
      for (int i = 0; i < degree_face; ++i) {
        face_descriptor fn = pmesh.face( CGAL::opposite(he, pmesh) );

        // if the neighbouring facet is regular
        if (check_regular_patch(fn, pmesh) == true) {
          // find vertices that determine the limit position of edge vertices
          std::vector<vertex_descriptor> vertices_edge_limit;
          halfedge_descriptor he_loop = he;
          for (int j = 0; j < degree_face; ++j) {
            vertex_descriptor vtx = CGAL::target(he_loop, pmesh);
            vertices_edge_limit.push_back(vtx);
            he_loop = CGAL::next(he_loop, pmesh);
          }

          typename PolygonMesh::size_type degree_fn = pmesh.degree(fn);
          he_loop = CGAL::opposite(he_loop, pmesh);
          he_loop = CGAL::next(he_loop, pmesh);
          for (int j = 0; j < degree_fn; ++j) {
            vertex_descriptor vtx = CGAL::target(he_loop, pmesh);
            if ( std::find(vertices_edge_limit.begin(), vertices_edge_limit.end(), vtx) == vertices_edge_limit.end() ) {
              vertices_edge_limit.push_back(vtx);
            }
            he_loop = CGAL::next(he_loop, pmesh);
          }

          // the limit coefficients
          double edge_limit_coeffs1[6] = {4./9., 1./9., 1./18., 2./9., 1./18., 1./9.};
          double edge_limit_coeffs2[6] = {2./9., 1./18., 1./9., 4./9., 1./9., 1./18.};

          std::vector<double> ep1(3, 0.), ep2(3, 0.);
          for (int j = 0; j < vertices_edge_limit.size(); ++j) {
            vertex_descriptor vtx = vertices_edge_limit[j];
            Kernel::Point_3 p = pmesh.point(vtx);

            for (int k = 0; k < 3; ++k) {
              ep1[k] += p[k]*edge_limit_coeffs1[j];
              ep2[k] += p[k]*edge_limit_coeffs2[j];
            }
          }

          Kernel::Point_3 ep_limit1(ep1[0], ep1[1], ep1[2]);
          Kernel::Point_3 ep_limit2(ep2[0], ep2[1], ep2[2]);

          if (i == 0) bezier_patch_ev[1] = ep_limit1, bezier_patch_ev[2] = ep_limit2;
          else if (i == 1) bezier_patch_ev[8] = ep_limit1, bezier_patch_ev[4] = ep_limit2;
          else if (i == 2) bezier_patch_ev[14] = ep_limit1, bezier_patch_ev[13] = ep_limit2;
          else if (i == 3) bezier_patch_ev[7] = ep_limit1, bezier_patch_ev[11] = ep_limit2;

        } // for edges with regular neighbouring facet
        else {
          std::vector< std::pair<int, int> > edge_vertices_LS;
          collect_edge_vertices_LS(edge_vertices_LS, i);

          // limit points on the edge to be approximated by least-square fitting
          std::vector< Kernel::Point_3 > edge_limit_points_LS;
          for (int j = 0; j < edge_vertices_LS.size(); ++j) {
            std::pair<int, int> edge_vertex_LS = edge_vertices_LS[j];
            Kernel::Point_3 limit_point = vtx_limit_map(edge_vertex_LS.first, edge_vertex_LS.second);
            edge_limit_points_LS.push_back(limit_point);
          }

          // least-square fitting by bezier points
          std::vector< Kernel::Point_3 > bezier_points;
          edge_least_square_fitting(bezier_points, edge_limit_points_LS);

          if (i == 0) bezier_patch_ev[2] = bezier_points[0], bezier_patch_ev[1] = bezier_points[1];
          else if (i == 1) bezier_patch_ev[4] = bezier_points[0], bezier_patch_ev[8] = bezier_points[1];
          else if (i == 2) bezier_patch_ev[13] = bezier_points[0], bezier_patch_ev[14] = bezier_points[1];
          else if (i == 3) bezier_patch_ev[11] = bezier_points[0], bezier_patch_ev[7] = bezier_points[1];

        } // for edges with irregular neighbouring facet

        he = CGAL::next(he, pmesh);

      } // loop over halfedges of the facet

      // centre bezier points within facet using least-square fitting
      std::vector< std::pair<int, int> > face_vertices_LS;
      collect_face_vertices_LS(face_vertices_LS);

      std::vector< Kernel::Point_3 > face_limit_points_LS;
      for (int j = 0; j < face_vertices_LS.size(); ++j) {
        std::pair<int, int> face_vertex_LS = face_vertices_LS[j];
        Kernel::Point_3 limit_point = vtx_limit_map(face_vertex_LS.first, face_vertex_LS.second);
        face_limit_points_LS.push_back(limit_point);
      }

      // least-square fitting by bezier points
      face_least_square_fitting(bezier_patch_ev, face_limit_points_LS);

      // rearrange Bezier control points
      std::map< std::pair<int, int>, Kernel::Point_3 > cps_bezier_ev;
      create_bezier_control_points(cps_bezier_ev, bezier_patch_ev);

      /*
            // output irregular bezier points
            for (int i = 0; i < bezier_patch_ev.size(); ++i) {
              Kernel::Point_3 cp = bezier_patch_ev[i];
              std::cout << cp << std::endl;
            }
       */

/*
      for (int j = 0; j < 11; ++j) {
        double v = 0. + 0.1*j;
        for (int i = 0; i < 11; ++i) {
          double u = 0. + 0.1*i;
          Kernel::Point_3 point = evaluate_bezier_surface(u, v, cps_bezier_ev);
          std::cout << point << std::endl;
        }
      }
*/
    } // irregular patches

  } // loop over facets in orginal mesh

/*
  Timer t;
  t.start();
  //Subdivision_method_3::CatmullClark_subdivision(pmesh, d);
  Subdivision_method_3::CatmullClark_subdivision(pmesh,
      Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, pmesh)).number_of_iterations(6));
  std::cerr << "Done (" << t.time() << " s)" << std::endl;

  std::ofstream out(out_file);
  out << pmesh;
*/
  return 0;
}

//=============================================================================
void onering_vertices_around_vertex(std::vector< typename boost::graph_traits<PolygonMesh>::vertex_descriptor >& onering_vertices_e,
                                    std::vector< typename boost::graph_traits<PolygonMesh>::vertex_descriptor >& onering_vertices_f,
                                    const typename boost::graph_traits<PolygonMesh>::vertex_descriptor& vtx,
                                    const PolygonMesh& pmesh)
{
  typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor he_vtx = pmesh.halfedge(vtx);

  // loop over one facet
  while (true) {
    halfedge_descriptor he_next = CGAL::next(he_vtx, pmesh);
    vertex_descriptor vtx_next = CGAL::target(he_next, pmesh);

    halfedge_descriptor he_next_peep = CGAL::next(he_next, pmesh);
    if ( CGAL::target(he_next_peep, pmesh) == vtx ) {
      onering_vertices_e.push_back(vtx_next);
      he_vtx = he_next_peep;
      break;
    }
    else {
      if ( CGAL::source(he_next, pmesh) == vtx ) onering_vertices_e.push_back(vtx_next);
      else onering_vertices_f.push_back(vtx_next);
      he_vtx = he_next;
    }
  }

  // loop over other incident facets
  bool is_continue = true;
  while (is_continue) {
    halfedge_descriptor he_opposite = CGAL::opposite(he_vtx, pmesh);
    while (true) {
      he_vtx = CGAL::next(he_opposite, pmesh);
      vertex_descriptor vtx_next = CGAL::target(he_vtx, pmesh);

      halfedge_descriptor he_next_peep = CGAL::next(he_vtx, pmesh);
      if ( CGAL::target(he_next_peep, pmesh) == vtx ) {
        std::vector<vertex_descriptor>::iterator it = std::find(onering_vertices_e.begin(), onering_vertices_e.end(), vtx_next);
        if (it != onering_vertices_e.end()) {
          is_continue = false;
          break;
        }
        else {
          onering_vertices_e.push_back(vtx_next);
          he_vtx = he_next_peep;
          break;
        }
      } // termination criteria
      else {
        if ( CGAL::source(he_vtx, pmesh) == vtx ) onering_vertices_e.push_back(vtx_next);
        else onering_vertices_f.push_back(vtx_next);
        he_opposite = he_vtx;
      }
    } // loop within one facet
  } // loop over facets

}

void collect_subdivision_map(std::map< typename boost::graph_traits<PolygonMesh>::face_descriptor, std::vector< typename boost::graph_traits<PolygonMesh>::halfedge_descriptor > >& subdiv_map,
                             const PolygonMesh& pmesh,
                             const PolygonMesh& pmesh_subdiv)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  typedef typename boost::graph_traits<PolygonMesh>::face_iterator face_iterator;

  face_iterator fitr = pmesh.faces_begin();
  for (fitr; fitr != pmesh.faces_end(); ++fitr) {
    face_descriptor f = *fitr;

    typename PolygonMesh::size_type degree_face = pmesh.degree(f);

    halfedge_descriptor he_min = pmesh.halfedge(f);
    halfedge_descriptor he_tmp = pmesh.halfedge(f);
    for (int i = 0; i < degree_face; ++i) {
      if (CGAL::target(he_tmp, pmesh) < CGAL::target(he_min, pmesh)) {
        he_min = he_tmp;
      }
      he_tmp = CGAL::next(he_tmp, pmesh);
    }

    halfedge_descriptor he = he_min;

    std::vector<halfedge_descriptor> halfedges_subdiv_face;

    for (int i = 0; i < degree_face; ++i) {
      halfedge_descriptor he_subdiv = he;

      // make sure halfedges pointing to targets in original mesh
      if ( CGAL::target(he_subdiv, pmesh_subdiv) != CGAL::target(he, pmesh) ) {
        he_subdiv = CGAL::next(he_subdiv, pmesh_subdiv);
        he_subdiv = CGAL::opposite(he_subdiv, pmesh_subdiv);
        he_subdiv = CGAL::next(he_subdiv, pmesh_subdiv);
      }

      halfedges_subdiv_face.push_back(he_subdiv);

      he = CGAL::next(he, pmesh);

    } // loop within original facet

    subdiv_map[f] = halfedges_subdiv_face;

  } // loop over original facets

}

void index_reorder_map(std::vector< std::pair<int, int> >& indices, const int& id1, const int& id2)
{
  std::vector<int> index1, index2;

  if (id1 == 0) {
    index1.push_back(0), index1.push_back(1), index1.push_back(2);
    index2.push_back(0), index2.push_back(1), index2.push_back(2);
  }
  else if (id1 == 1) {
    index1.push_back(0), index1.push_back(1), index1.push_back(2);
    index2.push_back(4), index2.push_back(3), index2.push_back(2);
  }
  else if (id1 == 2) {
    index1.push_back(4), index1.push_back(3), index1.push_back(2);
    index2.push_back(4), index2.push_back(3), index2.push_back(2);
  }
  else if (id1 == 3) {
    index1.push_back(4), index1.push_back(3), index1.push_back(2);
    index2.push_back(0), index2.push_back(1), index2.push_back(2);
  }

  std::vector< std::pair<int, int> > sub_indices;
  if (id2 == 0) {
    sub_indices.push_back( std::make_pair(0, 0) ), sub_indices.push_back( std::make_pair(1, 0) );
    sub_indices.push_back( std::make_pair(1, 1) ), sub_indices.push_back( std::make_pair(0, 1) );
  }
  else if (id2 == 1) {
    sub_indices.push_back( std::make_pair(2, 0) ), sub_indices.push_back( std::make_pair(2, 1) );
    sub_indices.push_back( std::make_pair(1, 1) ), sub_indices.push_back( std::make_pair(1, 0) );
  }
  else if (id2 == 2) {
    sub_indices.push_back( std::make_pair(2, 2) ), sub_indices.push_back( std::make_pair(1, 2) );
    sub_indices.push_back( std::make_pair(1, 1) ), sub_indices.push_back( std::make_pair(2, 1) );
  }
  else if (id2 == 3) {
    sub_indices.push_back( std::make_pair(0, 2) ), sub_indices.push_back( std::make_pair(0, 1) );
    sub_indices.push_back( std::make_pair(1, 1) ), sub_indices.push_back( std::make_pair(1, 2) );
  }

  for (int i = 0; i < sub_indices.size(); ++i) {
    std::pair<int, int> sub_index = sub_indices[i];
    int sub_id1 = index1[sub_index.first], sub_id2 = index2[sub_index.second];

    if (id1 == 0 || id1 == 2) indices.push_back( std::make_pair(sub_id1, sub_id2) );
    else if (id1 == 1 || id1 == 3) indices.push_back( std::make_pair(sub_id2, sub_id1) );

  }

}

void vertices_reorder_map(boost::numeric::ublas::bounded_matrix< vertex_descriptor, 5, 5 >& vtx_map,
                          const face_descriptor& f, const PolygonMesh& pmesh1, const PolygonMesh& pmesh2,
                          std::map< face_descriptor, std::vector< halfedge_descriptor > > subdiv_map_1,
                          std::map< face_descriptor, std::vector< halfedge_descriptor > > subdiv_map_2)
{
  // child facets of facet f
  std::vector<halfedge_descriptor> hes_subdiv_1 = subdiv_map_1[f];
  for (int i = 0; i < hes_subdiv_1.size(); ++i) {
    halfedge_descriptor he_subdiv_1 = hes_subdiv_1[i];
    face_descriptor f_subdiv_1 = pmesh1.face(he_subdiv_1);

    // child facets of facet f_subdiv_1
    std::vector<halfedge_descriptor> hes_subdiv_2 = subdiv_map_2[f_subdiv_1];
    for (int j = 0; j < hes_subdiv_2.size(); ++j) {
      halfedge_descriptor he_subdiv_2 = hes_subdiv_2[j];
      face_descriptor f_subdiv_2 = pmesh2.face(he_subdiv_2);

      std::vector< std::pair<int, int> > indices;
      index_reorder_map(indices, i, j);

      // vertices of facet f_subdiv_2 (reorder to start from the smallest index)
      typename PolygonMesh::size_type degree_f_subdiv_2 = pmesh2.degree(f_subdiv_2);

      halfedge_descriptor he_min = pmesh2.halfedge(f_subdiv_2);
      halfedge_descriptor he_tmp = he_min;

      for (int k = 0; k < degree_f_subdiv_2; ++k) {
        if ( CGAL::target(he_tmp, pmesh2) < CGAL::target(he_min, pmesh2) ) he_min = he_tmp;
        he_tmp = CGAL::next(he_tmp, pmesh2);
      }

      halfedge_descriptor he = he_min;

      for (int k = 0; k < degree_f_subdiv_2; ++k) {
        vertex_descriptor vtx = CGAL::target(he, pmesh2);
        std::pair<int, int> index = indices[k];
        vtx_map(index.first, index.second) = vtx;
        he = CGAL::next(he, pmesh2);
      }

    } // loop over child facets of a child facet
  } // loop over child facets of a facet in original mesh

}

void vertices_limit_map(boost::numeric::ublas::bounded_matrix< Kernel::Point_3, 5, 5 >& vtx_limit_map,
                        boost::numeric::ublas::bounded_matrix< vertex_descriptor, 5, 5 > vtx_map,
                        const PolygonMesh& pmesh2)
{
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      vertex_descriptor vtx = vtx_map(i, j);

      typename PolygonMesh::size_type degree_vtx = pmesh2.degree(vtx);

      // get one-ring neighbouring vertices of the vertex in pmesh2
      std::vector<vertex_descriptor> onering_vertices_e;
      std::vector<vertex_descriptor> onering_vertices_f;

      onering_vertices_around_vertex(onering_vertices_e, onering_vertices_f, vtx, pmesh2);

      std::vector<double> point(3, 0.);

      // centre vertex
      Kernel::Point_3 p_vtx = pmesh2.point(vtx);
      for (int ii = 0; ii < 3; ++ii) {
        point[ii] += p_vtx[ii]*degree_vtx/(degree_vtx + 5);
      }

      // edge vertices
      for (int ii = 0; ii < onering_vertices_e.size(); ++ii) {
        vertex_descriptor vtx_e = onering_vertices_e[ii];
        Kernel::Point_3 p_e = pmesh2.point(vtx_e);

        for (int jj = 0; jj < 3; ++jj) {
          point[jj] += p_e[jj]*4./(degree_vtx*(degree_vtx + 5));
        }
      }

      // vertices not incident to vtx
      for (int ii = 0; ii < onering_vertices_f.size(); ++ii) {
        vertex_descriptor vtx_f = onering_vertices_f[ii];
        Kernel::Point_3 p_f = pmesh2.point(vtx_f);

        for (int jj = 0; jj < 3; ++jj) {
          point[jj] += p_f[jj]/(degree_vtx*(degree_vtx + 5));
        }
      }

      Kernel::Point_3 point_limit(point[0], point[1], point[2]);

      vtx_limit_map(i, j) = point_limit;

    }
  } // loop over subdivided vertices of a facet
}

bool check_regular_patch(const face_descriptor& f, const PolygonMesh& pmesh)
{
  bool is_regular = true;

  typename PolygonMesh::size_type degree_face = pmesh.degree(f);

  halfedge_descriptor he = pmesh.halfedge(f);
  for (int i = 0; i < degree_face; ++i) {
    vertex_descriptor vtx = CGAL::target(he, pmesh);

    if (pmesh.degree(vtx) != 4) is_regular = false;

    he = CGAL::next(he, pmesh);
  }

  return is_regular;
}

void collect_edge_vertices_LS(std::vector< std::pair<int, int> >& indices, const int& edge_id)
{
  if (edge_id == 0) {
    indices.push_back( std::make_pair(0, 4) ), indices.push_back( std::make_pair(0, 3) );
    indices.push_back( std::make_pair(0, 2) ), indices.push_back( std::make_pair(0, 1) );
    indices.push_back( std::make_pair(0, 0) );
  }
  else if (edge_id == 1) {
    indices.push_back( std::make_pair(0, 0) ), indices.push_back( std::make_pair(1, 0) );
    indices.push_back( std::make_pair(2, 0) ), indices.push_back( std::make_pair(3, 0) );
    indices.push_back( std::make_pair(4, 0) );
  }
  else if (edge_id == 2) {
    indices.push_back( std::make_pair(4, 0) ), indices.push_back( std::make_pair(4, 1) );
    indices.push_back( std::make_pair(4, 2) ), indices.push_back( std::make_pair(4, 3) );
    indices.push_back( std::make_pair(4, 4) );
  }
  else if (edge_id == 3) {
    indices.push_back( std::make_pair(4, 4) ), indices.push_back( std::make_pair(3, 4) );
    indices.push_back( std::make_pair(2, 4) ), indices.push_back( std::make_pair(1, 4) );
    indices.push_back( std::make_pair(0, 4) );
  }
}

void edge_least_square_fitting(std::vector< Kernel::Point_3 >& bezier_points, const std::vector< Kernel::Point_3 >& edge_limit_points)
{
  Kernel::Vector_3 u(1./4., 1./2., 3./4.);

  Eigen::MatrixXd matrixL(3, 2);
  for (int i = 0; i < 3; ++i) {
    matrixL(i, 0) = 3.*u[i]*std::pow(1. - u[i], 2.);
    matrixL(i, 1) = 3.*std::pow(u[i], 2.)*(1. - u[i]);
  }

  Eigen::MatrixXd matrixLT = matrixL.transpose();

  Eigen::MatrixXd matrixP(3, 2);
  for (int i = 0; i < 3; ++i) {
    matrixP(i, 0) = std::pow(1. - u[i], 3.);
    matrixP(i, 1) = std::pow(u[i], 3.);
  }

  Eigen::MatrixXd matrixB1(3, 3);
  for (int i = 0; i < 3; ++i) {
    Kernel::Point_3 p = edge_limit_points[i + 1];
    for (int j = 0; j < 3; ++j) {
      matrixB1(i, j) = p[j];
    }
  }

  Eigen::MatrixXd matrixB2(2, 3);
  for (int i = 0; i < 3; ++i) {
    Kernel::Point_3 p1 = edge_limit_points[0];
    Kernel::Point_3 p2 = edge_limit_points[4];

    matrixB2(0, i) = p1[i];
    matrixB2(1, i) = p2[i];
  }

  Eigen::MatrixXd matrixR = matrixB1 - matrixP*matrixB2;

  Eigen::MatrixXd matrix_tmp = matrixLT*matrixL;
  Eigen::MatrixXd matrix_tmp_inv = matrix_tmp.inverse();
  Eigen::MatrixXd matrixLT_tmp = matrix_tmp_inv*matrixLT;

  Eigen::MatrixXd result = matrixLT_tmp*matrixR;

  Kernel::Point_3 bezier_point1( result(0, 0), result(0, 1), result(0, 2) );
  Kernel::Point_3 bezier_point2( result(1, 0), result(1, 1), result(1, 2) );

  bezier_points.push_back(bezier_point1);
  bezier_points.push_back(bezier_point2);

}

void collect_face_vertices_LS(std::vector< std::pair<int, int> >& indices)
{
  indices.push_back( std::make_pair(1, 1) ), indices.push_back( std::make_pair(1, 2) ), indices.push_back( std::make_pair(1, 3) );
  indices.push_back( std::make_pair(2, 1) ), indices.push_back( std::make_pair(2, 2) ), indices.push_back( std::make_pair(2, 3) );
  indices.push_back( std::make_pair(3, 1) ), indices.push_back( std::make_pair(3, 2) ), indices.push_back( std::make_pair(3, 3) );
}

void face_least_square_fitting(std::map<int, Kernel::Point_3>& bezier_points, const std::vector< Kernel::Point_3 >& face_limit_points)
{
  std::vector< std::pair<double, double> > params(9);
  params[0] = std::make_pair(1./4., 1./4.), params[1] = std::make_pair(1./2., 1./4.), params[2] = std::make_pair(3./4., 1./4.);
  params[3] = std::make_pair(1./4., 1./2.), params[4] = std::make_pair(1./2., 1./2.), params[5] = std::make_pair(3./4., 1./2.);
  params[6] = std::make_pair(1./4., 3./4.), params[7] = std::make_pair(1./2., 3./4.), params[8] = std::make_pair(3./4., 3./4.);

  Eigen::MatrixXd bernstein_full(9, 16); // 16 basis functions evaluated at 9 parameters
  for (int i = 0; i < 9; ++i) {
    double u = params[i].first, v = params[i].second;

    Eigen::Vector4d u_bernstein, v_bernstein;
    u_bernstein(0) = std::pow(1. - u, 3.), u_bernstein(1) = 3.*u*std::pow(1. - u, 2.), u_bernstein(2) = 3.*u*u*(1. - u), u_bernstein(3) = std::pow(u, 3.);
    v_bernstein(0) = std::pow(1. - v, 3.), v_bernstein(1) = 3.*v*std::pow(1. - v, 2.), v_bernstein(2) = 3.*v*v*(1. - v), v_bernstein(3) = std::pow(v, 3.);

    int count = 0;
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        bernstein_full(i, count) = u_bernstein(k)*v_bernstein(j);
        count += 1;
      }
    }
  }

  int labelL[4] = {5, 6, 9, 10};
  int labelR[12] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};

  Eigen::MatrixXd bernsteinL(9, 4), bernsteinR(9, 12);
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 4; ++j) {
      int label = labelL[j];
      bernsteinL(i, j) = bernstein_full(i, label);
    }

    for (int j = 0; j < 12; ++j) {
      int label = labelR[j];
      bernsteinR(i, j) = bernstein_full(i, label);
    }
  }

  Eigen::MatrixXd bezierR(12, 3);
  for (int i = 0; i < 12; ++i) {
    int label = labelR[i];

    for (int j = 0; j < 3; ++j) {
      bezierR(i, j) = bezier_points[label][j];
    }
  }

  Eigen::MatrixXd targets(9, 3);
  for (int i = 0; i < 9; ++i) {
    Kernel::Point_3 point = face_limit_points[i];
    for (int j = 0; j < 3; ++j) {
      targets(i, j) = point[j];
    }
  }

  Eigen::MatrixXd matrixR = targets - bernsteinR*bezierR;

  Eigen::MatrixXd bernsteinLT = bernsteinL.transpose();
  Eigen::MatrixXd matrix_tmp = bernsteinLT*bernsteinL;
  Eigen::MatrixXd matrix_tmp_inv = matrix_tmp.inverse();
  Eigen::MatrixXd matrixLT_tmp = matrix_tmp_inv*bernsteinLT;

  Eigen::MatrixXd result = matrixLT_tmp*matrixR;

  for (int i = 0; i < 4; ++i) {
    int label = labelL[i];

    Kernel::Point_3 point( result(i, 0), result(i, 1), result(i, 2) );
    bezier_points[label] = point;
  }

}

void bspline_to_bezier(std::map<int, Kernel::Point_3>& bezier_points, const Eigen::MatrixXd& cps_spline)
{
  const int degree = 3;
  const int num_segments = degree + 1;

  Eigen::MatrixXd coeffs_bezier;
  coeffs_bezier.resize(num_segments, degree + 1);

  // Bezier coefficients for a cubic b-spline
  coeffs_bezier(0, 0) = 1., coeffs_bezier(0, 1) = 0., coeffs_bezier(0, 2) = 0., coeffs_bezier(0, 3) = 0.;
  coeffs_bezier(1, 0) = 4., coeffs_bezier(1, 1) = 4., coeffs_bezier(1, 2) = 2., coeffs_bezier(1, 3) = 1.;
  coeffs_bezier(2, 0) = 1., coeffs_bezier(2, 1) = 2., coeffs_bezier(2, 2) = 4., coeffs_bezier(2, 3) = 4.;
  coeffs_bezier(3, 0) = 0., coeffs_bezier(3, 1) = 0., coeffs_bezier(3, 2) = 0., coeffs_bezier(3, 3) = 1.;

  coeffs_bezier /= 6.;

  // transformation matrix from B-spline to Bezier
  Eigen::MatrixXd coeffs_transform;
  coeffs_transform.resize( (degree + 1)*num_segments, num_segments*num_segments );

  for (int j = 0; j < num_segments; ++j) {
    Eigen::MatrixXd coeffsj = coeffs_bezier.block(j, 0, 1, degree + 1);

    for (int i = 0; i < num_segments; ++i) {
      Eigen::MatrixXd coeffsi = coeffs_bezier.block(i, 0, 1, degree + 1);

      int col_id = i + j*num_segments;

      for (int jj = 0; jj < degree + 1; ++jj) {
        for (int ii = 0; ii < degree + 1; ++ii) {
          double coeff = coeffsi(0, ii)*coeffsj(0, jj);
          int row_id = ii + jj*(degree + 1);
          coeffs_transform(row_id, col_id) = coeff;
        }
      }
    }
  }

  // Bezier control points
  Eigen::MatrixXd cps_bezier = coeffs_transform*cps_spline;

  for (int i = 0; i < num_segments*num_segments; ++i) {
    Kernel::Point_3 point( cps_bezier(i, 0), cps_bezier(i, 1), cps_bezier(i, 2) );
    bezier_points[i] = point;
  }
}

void onering_vertices_around_face(std::vector< vertex_descriptor >& onering_vertices,
                                  const face_descriptor& f,
                                  const PolygonMesh& pmesh)
{
  std::vector< vertex_descriptor > vertices_face;
  std::vector< vertex_descriptor > neighbour_vertices_face;

  typename PolygonMesh::size_type degree_face = pmesh.degree(f);
  halfedge_descriptor he = pmesh.halfedge(f);
  for (int i = 0; i < degree_face; ++i) {
    vertex_descriptor vtx = CGAL::target(he, pmesh);
    vertices_face.push_back(vtx);
    neighbour_vertices_face.push_back(vtx);
    he = CGAL::next(he, pmesh);
  }

  bool is_continue = true;
  while (is_continue) {
    halfedge_descriptor he_opposite = CGAL::opposite(he, pmesh);
    while (true) {
      he = CGAL::next(he_opposite, pmesh);
      vertex_descriptor vtx = CGAL::target(he, pmesh);

      std::vector<vertex_descriptor>::iterator it = std::find(vertices_face.begin(), vertices_face.end(), vtx);

      if (it != vertices_face.end()) {
        break;
      }
      else {
        neighbour_vertices_face.push_back(vtx);
        he_opposite = he;

        if (neighbour_vertices_face.size() == 16) {
          is_continue = false;
          break;
        }
      }
    }
  }

  // define an index array
  int id_array[16] = {5, 9, 10, 6, 2, 1, 0, 4, 8, 12, 13, 14, 15, 11, 7, 3};

  for (int i = 0; i < neighbour_vertices_face.size(); ++i) {
    int id = id_array[i];
    onering_vertices[id] = neighbour_vertices_face[i];
  }

}

void create_bezier_control_points(std::map< std::pair<int, int>, Kernel::Point_3 >& cps_bezier, std::map<int, Kernel::Point_3> bezier_patch)
{
  int d1 = 3, d2 = 3;

  int count = 0;

  // note the ordering of control points (related to u and v directions)
  for (int j = 0; j < d2 + 1; ++j) {
    for (int i = 0; i < d1 + 1; ++i) {
      std::pair<int, int> id = std::make_pair(i, j);
      Kernel::Point_3 point = bezier_patch[count];
      cps_bezier[id] = point;
      count += 1;
    }
  }

}

Kernel::Point_3 evaluate_bezier_surface(const double& u, const double& v,
                                        std::map< std::pair<int, int>, Kernel::Point_3 > cps_bezier)
{
  std::vector<double> p(3, 0.);

  int d1 = 3, d2 = 3;

  for (int j = 0; j < d2 + 1; ++j) {
    double coeff1 = boost::math::binomial_coefficient<double>(d2, j);
    double bernstein1 = coeff1*std::pow(v, j)*std::pow(1. - v, d2 - j);

    for (int i = 0; i < d1 + 1; ++i) {
      double coeff2 = boost::math::binomial_coefficient<double>(d1, i);
      double bernstein2 = coeff2*std::pow(u, i)*std::pow(1. - u, d1 -i);

      std::pair<int, int> id = std::make_pair(i, j);
      Kernel::Point_3 bezier_point = cps_bezier[id];

      for (int k = 0; k < 3; ++k) {
        p[k] += bernstein1*bernstein2*bezier_point[k];
      }
    }

  }

  Kernel::Point_3 point( p[0], p[1], p[2] );

  return point;
}
