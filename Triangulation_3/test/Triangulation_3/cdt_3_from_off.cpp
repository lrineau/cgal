//#define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#define NO_TRY_CATCH 1
// #define CGAL_DEBUG_CDT_3 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <ranges>
#include <optional>

#if NO_TRY_CATCH
# define CDT_3_try      if (true)
# define CDT_3_catch(X) if (false)
# define CDT_3_throw_exception_again
#else
// Else proceed normally.
# define CDT_3_try      try
# define CDT_3_catch(X) catch(X)
# define CDT_3_throw_exception_again throw
#endif

#if CGAL_CDT_3_USE_EPECK

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;

#else // use Epick

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

#endif // use Epick

using Vb = CGAL::Base_with_time_stamp<CGAL::Constrained_Delaunay_triangulation_vertex_base_3<K>>;
using Cb = CGAL::Constrained_Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
using CDT = CGAL::Constrained_Delaunay_triangulation_3<Delaunay>;
using Point = Delaunay::Point;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor   = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor   = boost::graph_traits<Mesh>::face_descriptor;

struct CDT_options
{
  bool merge_facets = false;
  double ratio = 0.;
  std::string input_filename = CGAL::data_file_path("meshes/mpi.off");
  std::string output_filename{"dump.off"};
  std::string dump_patches_after_merge_filename{};
  std::string dump_patches_borders_prefix{};
  std::string dump_after_conforming_filename{};
};

int go(Mesh, CDT_options);

void help(std::ostream& out) {
  out << R"(
Usage: cdt_3_from_off [options] input.off output.off

  input.off: input mesh
  output.off: output mesh

  --merge-facets: merge facets into patches
  --ratio: ratio of faces to remove (default: 0)
  --dump-patches-after-merge: dump patches after merging facets
  --dump-patches-borders-prefix: dump patches borders
  --dump-after-conforming: dump mesh after conforming
)";
}

int main(int argc, char* argv[])
{
  std::cerr.precision(17);
  std::cout.precision(17);

  CDT_options options;
  int positional = 0;

  for(int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if(arg == "--merge-facets") {
      options.merge_facets = true;
    } else if(arg == "--ratio") {
      assert(i + 1 < argc);
      options.ratio = std::stod(argv[++i]);
    } else if(arg == "--dump-patches-after-merge") {
      assert(i + 1 < argc);
      options.dump_patches_after_merge_filename = argv[++i];
    } else if(arg == "--dump-patches-borders-prefix") {
      assert(i + 1 < argc);
      options.dump_patches_borders_prefix = argv[++i];
    } else if(arg == "--dump-after-conforming") {
      assert(i + 1 < argc);
      options.dump_after_conforming_filename = argv[++i];
    } else if(arg == "--help") {
      help(std::cout);
      return 0;
    } else if(arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << '\n';
      help(std::cerr);
      return 1;
    } else {
      switch(positional) {
        case 0:
          options.input_filename = arg;
          ++positional;
          break;
        case 1:
          options.output_filename = arg;
          ++positional;
          break;
        case 2:
          options.ratio = std::stod(arg);
          ++positional;
          break;
        default:
          std::cerr << "Too many arguments\n";
          return 1;
      }
    }
  }

  Mesh mesh;
  const bool ok = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(options.input_filename, mesh);
  if (!ok)
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  if(options.ratio == 0.) {
    return go(std::move(mesh), std::move(options));
  }
  auto nb_buckets = static_cast<int>(std::floor(1 / options.ratio)) + 1;
  std::cerr << "RATIO: " << options.ratio << '\n';

  const Mesh orig_mesh{mesh};
  Mesh bad_mesh{mesh};
  for(int bucket = 0; bucket < nb_buckets;) {
    const auto nb_faces = mesh.number_of_faces();
    auto nb_to_skip = static_cast<int>(std::round(nb_faces * options.ratio));
    if(nb_to_skip < 1) {
      nb_to_skip = 1;
      nb_buckets = nb_faces;
    }
    if(bucket == 0) {
      std::cerr << "NB BUCKETS: " << nb_buckets << '\n';
    }

    auto simplify = [&](Mesh& m) {
      std::cerr << "nb_to_skip: " << nb_to_skip << '\n';
      std::cerr << "bucket: " << bucket << '\n';
      const auto start = (std::min)(bucket * nb_to_skip, static_cast<int>(m.number_of_faces()));
      const auto end = (std::min)(start + nb_to_skip, static_cast<int>(m.number_of_faces()));
      std::cerr << "SKIP from " << start << " to " << end << '\n';
      for(auto i = end - 1; i >= start; --i) {
        const auto f = m.faces().begin() + i;
        CGAL::Euler::remove_face(halfedge(*f, m), m);
      }
      assert(m.is_valid(true));
      std::cerr << "number of faces: " << m.number_of_faces() << '\n';
      if(m.number_of_faces() >= nb_faces) {
        std::cerr << "ERROR: could not simplify mesh\n";
        std::exit(EXIT_FAILURE);
      }
    };

    simplify(mesh);
        std::ofstream current("current_mesh.off");
        current.precision(17);
        current << mesh;
        current.close();

    try {
      go(mesh, options);
    } catch(CGAL::Failure_exception& e) {
      if(e.expression().find(std::string("orientation")) != std::string::npos)
      {
        std::cerr << "BAD MESH! " << mesh.number_of_faces() << " faces\n";
        std::ofstream bad("bad_mesh.off");
        bad.precision(17);
        bad << mesh;
        bad_mesh = mesh;
        bucket = 0;
        continue;
      } else {
        std::cerr << "ERROR MESH: " << e.what() << '\n';
        std::ofstream error("error_mesh.off");
        error.precision(17);
        error << mesh;
        std::cerr << "go on...\n";
      }
    }
    std::cerr << "GOOD MESH :-( " << mesh.number_of_faces() << " faces\n";
    mesh = bad_mesh;
    ++bucket;
  }
}

template <typename Range_of_segments>
auto segment_soup_to_polylines(Range_of_segments&& segment_soup) {
  using Point = decltype([&](){ using std::begin; auto [a, b] = *begin(segment_soup); return a; } ());

  std::vector<std::vector<Point>> polylines;

  using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Point>;
  using Map2v = std::map<Point, typename Graph::vertex_descriptor>;
  Graph graph;
  Map2v map2v;
  auto get_v = [&](const Point& p) {
    auto it = map2v.find(p);
    if(it != map2v.end()) return it->second;
    auto v = boost::add_vertex(p, graph);
    map2v.emplace(p, v);
    return v;
  };
  for(auto [a, b]: segment_soup) {
    auto va = get_v(a);
    auto vb = get_v(b);
    boost::add_edge(va, vb, graph);
  }

  struct Polylines_visitor
  {
    Graph& graph;
    std::vector<std::vector<Point>>& polylines;

    void start_new_polyline() { polylines.emplace_back(); }
    void add_node(typename Graph::vertex_descriptor vd) { polylines.back().push_back(graph[vd]); }
    void end_polyline() {}
  };
  Polylines_visitor visitor{graph, polylines};
  CGAL::split_graph_into_polylines(graph, visitor);

  return polylines;
}

int go(Mesh mesh, CDT_options options) {
  CDT cdt;
  auto pmap = get(CGAL::vertex_point, mesh);

  auto [patch_id_map, ok] = mesh.add_property_map<face_descriptor, int>("f:patch_id", -1);
  assert(ok); CGAL_USE(ok);
  auto [v_selected_map, ok2] = mesh.add_property_map<vertex_descriptor, bool>("v:selected", false);
  assert(ok2); CGAL_USE(ok2);
  int nb_patches = 0;
  std::vector<std::vector<std::pair<vertex_descriptor, vertex_descriptor>>> patch_edges;
  if(options.merge_facets) {
    for(auto f: faces(mesh))
    {
      if(get(patch_id_map, f) >= 0) continue;
      patch_edges.emplace_back();
      auto& edges = patch_edges.back();
      std::stack<face_descriptor> f_stack;
      f_stack.push(f);
      while(!f_stack.empty()) {
        auto f = f_stack.top();
        f_stack.pop();
        if(get(patch_id_map, f) >= 0) continue;
        put(patch_id_map, f, nb_patches);
        for(auto h: CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
          auto opp = opposite(h, mesh);
          if(is_border_edge(opp, mesh)) {
            auto va = source(h, mesh);
            auto vb = target(h, mesh);
            edges.emplace_back(va, vb);
            put(v_selected_map, va, true);
            put(v_selected_map, vb, true);
            continue;
          }
          auto n = face(opp, mesh);
          auto a = get(pmap, source(h, mesh));
          auto b = get(pmap, target(h, mesh));
          auto c = get(pmap, target(next(h, mesh), mesh));
          auto d = get(pmap, target(next(opp, mesh), mesh));
          if(CGAL::orientation(a, b, c, d) != CGAL::COPLANAR) {
            auto va = source(h, mesh);
            auto vb = target(h, mesh);
            edges.emplace_back(va, vb);
            put(v_selected_map, va, true);
            put(v_selected_map, vb, true);
            continue;
          }
          if(get(patch_id_map, n) >= 0) continue;
          f_stack.push(n);
        }
      }
      ++nb_patches;
    }
    if(!options.dump_patches_after_merge_filename.empty()) {
      std::ofstream out(options.dump_patches_after_merge_filename);
      CGAL::IO::write_PLY(out, mesh);
    }
  }
  if(!options.dump_patches_borders_prefix.empty()) {
    for(int i = 0; i < nb_patches; ++i) {
      std::stringstream ss;
      ss << options.dump_patches_borders_prefix << i << ".polylines.txt";
      std::ofstream out(ss.str());
      out.precision(17);
      const auto& edges = patch_edges[i];
      std::cerr << "Patch p#" << i << " has " << edges.size() << " edges\n";
      const auto polylines = segment_soup_to_polylines(edges);
      for(const auto& polyline: polylines) {
        out << polyline.size() << "    ";
        for(auto v: polyline) {
          out << get(pmap, v) << "  ";
        }
        out << '\n';
      }
      out.close();
      std::cerr << "  " << polylines.size() << " polylines\n";
      for(const auto& polyline: polylines) {
        std::cerr << "    - " << polyline.size() << " vertices\n";
        assert(polyline.front() == polyline.back());
      }
    }
  }

  int exit_code = EXIT_SUCCESS;

  auto finally = [&cdt, &options]() {
    {
      std::ofstream dump("dump.binary.cgal");
      CGAL::IO::save_binary_file(dump, cdt);
    }
    {
      std::ofstream dump(options.output_filename);
      dump.precision(17);
      cdt.write_facets(dump, cdt, std::views::filter(cdt.finite_facets(), [&](auto f) {
          return cdt.is_constrained(f);
      }));
    }
    {
      std::ofstream missing_faces("dump_missing_faces.polylines.txt");
      missing_faces.precision(17);
      cdt.recheck_constrained_Delaunay();
      if(cdt.write_missing_subfaces_file(missing_faces)) {
        std::cerr << "ERROR: Missing subfaces!\n";
      }
    }
    {
      std::ofstream missing_edges("dump_missing_segments.polylines.txt");
      missing_edges.precision(17);
      if(cdt.write_missing_segments_file(missing_edges)) {
        std::cerr << "ERROR: Missing segments!\n";
      }
    }
  };

  for(auto v: vertices(mesh)) {
    if(options.merge_facets && false == get(v_selected_map, v)) continue;
    cdt.insert(get(pmap, v));
  }
  if(cdt.dimension() < 3) {
    const auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    double d_x = bbox.xmax() - bbox.xmin();
    double d_y = bbox.ymax() - bbox.ymin();
    double d_z = bbox.zmax() - bbox.zmin();

    const double max_d = (std::max)(d_x, (std::max)(d_y, d_z));
    if(d_x == 0) d_x = max_d;
    if(d_y == 0) d_y = max_d;
    if(d_z == 0) d_z = max_d;

    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));
  }
  int poly_id = 0;
  CDT_3_try {
    if(options.merge_facets) {
      for(int i = 0; i < nb_patches; ++i) {
        auto& edges = patch_edges[i];
        auto polylines = segment_soup_to_polylines(edges);
        while(true) {
          const auto non_closed_polylines_begin =
              std::partition(polylines.begin(), polylines.end(),
                             [](const auto& polyline) { return polyline.front() == polyline.back(); });
          if(non_closed_polylines_begin == polylines.end())
            break;
          edges.clear();
          for(auto it = non_closed_polylines_begin; it != polylines.end(); ++it) {
            auto& polyline = *it;
            for(auto it = polyline.begin(), end = polyline.end() - 1; it != end; ++it) {
              edges.emplace_back(*it, *(it + 1));
            }
          }
          polylines.erase(non_closed_polylines_begin, polylines.end());
          auto other_polylines = segment_soup_to_polylines(edges);
          polylines.insert(polylines.end(),
                           std::make_move_iterator(other_polylines.begin()),
                           std::make_move_iterator(other_polylines.end()));
        }

        std::optional<int> face_index;
        for(auto& polyline: polylines) {
          assert(polyline.front() == polyline.back());
          polyline.pop_back();
          if(face_index) {
            cdt.insert_constrained_polygon(
              polyline | std::views::transform([&](vertex_descriptor v) { return get(pmap, v); }),
              false,
              *face_index);
          } else {
            face_index = cdt.insert_constrained_polygon(
              polyline | std::views::transform([&](vertex_descriptor v) { return get(pmap, v); }),
              false);
          }
        }
      }
    } else {
      for(auto face_descriptor : faces(mesh)) {
        std::vector<Point_3> polygon;
        const auto he = halfedge(face_descriptor, mesh);
        for(auto vertex_it : CGAL::vertices_around_face(he, mesh)) {
          polygon.push_back(get(pmap, vertex_it));
        }
  #if CGAL_DEBUG_CDT_3
        std::cerr << "NEW POLYGON #" << poly_id << '\n';
  #endif // CGAL_DEBUG_CDT_3
        const auto coplanar = polygon.size() < 3 ||
            std::all_of(polygon.begin(), polygon.end(),
                        [p1 = polygon[0], p2 = polygon[1], p3 = polygon[2]](auto p) {
                          const auto coplanar =
                              CGAL::orientation(p1, p2, p3, p) == CGAL::COPLANAR;
                          if(!coplanar) {
                            std::cerr << "Non coplanar points: " << p1 << ", " << p2
                                      << ", " << p3 << ", " << p << '\n'
                                      << "  volume: " << volume(p1, p2, p3, p) << '\n';

                          }
                          return coplanar;
                        });
        if(!coplanar) {
          std::ofstream out(std::string("dump_noncoplanar_polygon_") + std::to_string(poly_id) + ".off");
          out.precision(17);
          out << "OFF\n" << polygon.size() << " 1 0\n";
          for(auto p : polygon) {
            out << p << '\n';
          }
          out << polygon.size() << ' ';
          for(std::size_t i = 0u, end = polygon.size(); i < end; ++i) {
            out << ' ' << i;
          }
          out << '\n';
          std::cerr << "Polygon is not coplanar\n";
        }
        try {
          [[maybe_unused]] auto id = cdt.insert_constrained_polygon(polygon, false);
          assert(id == poly_id);
          ++poly_id;
        } catch(int error) {
          exit_code = error;
        }
        // std::ofstream dump("dump.binary.cgal");
        // CGAL::Mesh_3::save_binary_file(dump, cdt);
      }
    } // not merge_facets
    cdt.restore_Delaunay();

    if(!options.dump_after_conforming_filename.empty()) {
      for(auto e: edges(mesh)) {
        auto he = halfedge(e, mesh);
        auto vd1 = target(he, mesh);
        auto vd2 = source(he, mesh);
        if(!get(v_selected_map, vd1) || !get(v_selected_map, vd2)) continue;
        auto p1 = get(pmap, vd1);
        auto p2 = get(pmap, vd2);
        auto n = cdt.number_of_vertices();
        auto v1 = cdt.insert(p1);
        auto v2 = cdt.insert(p2);
        std::cerr << std::format("edge  {}  {}\n", CGAL::IO::oformat(p1), CGAL::IO::oformat(p2));
        CGAL_assertion(n == cdt.number_of_vertices());
        auto steiner_vertices = cdt.sequence_of_Steiner_vertices(v1, v2);
        if(!steiner_vertices) continue;
        for(auto v: *steiner_vertices) {
          he = CGAL::Euler::split_edge(he, mesh);
          put(pmap, target(he, mesh), v->point());
        }
      }
      std::ofstream out_mesh(options.dump_after_conforming_filename);
      out_mesh.precision(17);
      out_mesh << mesh;
      out_mesh.close();
    }

    std::cerr << "Number of vertices after conforming: " << cdt.number_of_vertices() << '\n';
    assert(cdt.is_conforming());
    assert(cdt.Delaunay::is_valid(true));
    assert(cdt.is_valid(true));
    if(exit_code == EXIT_SUCCESS) {
      try {
        cdt.restore_constrained_Delaunay();
        std::cerr << "Number of vertices after CDT: " << cdt.number_of_vertices() << '\n';
      } catch(int error) {
        exit_code = error;
      }
    }
  } CDT_3_catch(CGAL::Failure_exception&) {
    finally();
    CDT_3_throw_exception_again;
  }
  finally();
  assert(cdt.is_conforming());
  assert(cdt.is_valid(true));


  return exit_code;
}