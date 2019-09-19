// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_INTERNAL_TET_REMESHING_HELPERS_H
#define CGAL_INTERNAL_TET_REMESHING_HELPERS_H

#include <utility>

#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/triangulation_3_helpers.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  enum Subdomain_relation { EQUAL, DIFFERENT, INCLUDED, INCLUDES };
  enum Sliver_removal_result { INVALID_ORIENTATION, INVALID_CELL, INVALID_VERTEX,
    NOT_FLIPPABLE, EDGE_PROBLEM, VALID_FLIP, NO_BEST_CONFIGURATION, EXISTING_EDGE };


  namespace helpers
  {
    template<typename C3T3, typename CellSelector>
    bool is_boundary(const C3T3& c3t3,
      const typename C3T3::Triangulation::Edge& e,
      CellSelector cell_selector)
    {
      typedef typename C3T3::Triangulation   Tr;
      typedef typename Tr::Facet_circulator  Facet_circulator;
      typedef typename Tr::Facet             Facet;

      Facet_circulator fcirc = c3t3.triangulation().incident_facets(e);
      Facet_circulator fend = fcirc;
      std::vector<Facet> boundary_facets;

      do
      {
        Facet f = *fcirc;
        if (c3t3.is_in_complex(f))
          return true;
        else if (cell_selector(f.first) // XOR
          ^ cell_selector(f.first->neighbor(f.second)))
          return true;
        else if (c3t3.triangulation().is_infinite(f) //XOR
          ^ c3t3.triangulation().is_infinite(f.first->neighbor(f.second)))
          return true;

        ++fcirc;
      } while (fcirc != fend);

      return false;
    }

    template<typename C3t3, typename CellSelector>
    bool is_boundary_edge(const typename C3t3::Vertex_handle& v0,
      const typename C3t3::Vertex_handle& v1,
      const C3t3& c3t3,
      CellSelector cell_selector)
    {
      typedef typename C3t3::Edge        Edge;
      typedef typename C3t3::Cell_handle Cell_handle;

      Cell_handle cell;
      int i0, i1;
      if (c3t3.triangulation().tds().is_edge(v0, v1, cell, i0, i1))
        return is_boundary(c3t3, Edge(cell, i0, i1), cell_selector);
      else
        return false;
    }

    template<typename C3t3, typename CellSelector>
    bool is_boundary_vertex(const typename C3t3::Vertex_handle& v,
      const C3t3& c3t3,
      CellSelector cell_selector)
    {
      typedef typename C3t3::Facet Facet;
      std::vector<Facet> facets;
      c3t3.triangulation().incident_facets(v, std::back_inserter(facets));

      BOOST_FOREACH(Facet f, facets)
      {
        if (c3t3.is_in_complex(f))
          return true;
        if (cell_selector(f.first) ^ cell_selector(f.first->neighbor(f.second)))
          return true;
      }
      return false;
    }

    template<typename C3t3, typename CellSelector>
    bool is_edge_in_complex(const typename C3t3::Vertex_handle& v0,
      const typename C3t3::Vertex_handle& v1,
      const C3t3& c3t3,
      CellSelector /*cell_selector*/)
    {
      typedef typename C3t3::Edge        Edge;
      typedef typename C3t3::Cell_handle Cell_handle;

      Cell_handle cell;
      int i0, i1;
      if (c3t3.triangulation().tds().is_edge(v0, v1, cell, i0, i1))
        return c3t3.is_in_complex(Edge(cell, i0, i1));
      else
        return false;
    }

    template<typename C3t3, typename CellSelector>
    bool topology_test(const typename C3t3::Edge& edge,
      const C3t3& c3t3,
      CellSelector cell_selector)
    {
      typedef typename C3t3::Vertex_handle Vertex_handle;
      typedef typename C3t3::Triangulation::Facet_circulator Facet_circulator;
      typedef typename C3t3::Subdomain_index Subdomain_index;

      Vertex_handle v0 = edge.first->vertex(edge.second);
      Vertex_handle v1 = edge.first->vertex(edge.third);

      Facet_circulator fcirc = c3t3.triangulation().incident_facets(edge);
      Facet_circulator fdone = fcirc;
      do
      {
        if (c3t3.triangulation().is_infinite(fcirc->first))
          continue;

        Subdomain_index si_circ = fcirc->first->subdomain_index();
        Subdomain_index si_neigh = fcirc->first->neighbor(fcirc->second)->subdomain_index();
        if (si_circ == si_neigh)
        {
          //Get the ids of the opposite vertices
          for (int i = 1; i < 4; i++)
          {
            Vertex_handle vi = fcirc->first->vertex((fcirc->second + i) % 4);
            if (vi != v0 && vi != v1 && nb_incident_subdomains(vi, c3t3) > 1)
            {
              if (is_edge_in_complex(v0, vi, c3t3, cell_selector)
                && is_edge_in_complex(v1, vi, c3t3, cell_selector))
                return false;
            }
          }
        }
      } while (++fcirc != fdone);

      return true;
    }

    template<typename C3t3>
    Subdomain_relation compare_subdomains(typename C3t3::Vertex_handle v0,
      typename C3t3::Vertex_handle v1,
      const C3t3& c3t3)
    {
      typedef typename C3t3::Subdomain_index Subdomain_index;

      std::vector<Subdomain_index> subdomains_v0;
      incident_subdomains(v0, c3t3, std::back_inserter(subdomains_v0));
      std::sort(subdomains_v0.begin(), subdomains_v0.end());

      std::vector<Subdomain_index> subdomains_v1;
      incident_subdomains(v1, c3t3, std::back_inserter(subdomains_v1));
      std::sort(subdomains_v1.begin(), subdomains_v1.end());

      if (subdomains_v0.size() == subdomains_v1.size())
      {
        for (unsigned int i = 0; i < subdomains_v0.size(); i++)
          if (subdomains_v0[i] != subdomains_v1[i])
            return DIFFERENT;
        return EQUAL;
      }
      else
      {
        std::vector<Subdomain_index>
          intersection((std::min)(subdomains_v0.size(), subdomains_v1.size()), -1);
        typename std::vector<Subdomain_index>::iterator
          end_it = std::set_intersection(subdomains_v0.begin(), subdomains_v0.end(),
            subdomains_v1.begin(), subdomains_v1.end(),
            intersection.begin());
        std::ptrdiff_t intersection_size = (end_it - intersection.begin());

        if (subdomains_v0.size() > subdomains_v1.size()
          && intersection_size == std::ptrdiff_t(subdomains_v1.size()))
        {
          return INCLUDES;
        }
        else if (intersection_size == std::ptrdiff_t(subdomains_v0.size())) {
          return INCLUDED;
        }
      }
      return DIFFERENT;
    }



    template<typename C3t3, typename CellSelector>
    void get_edge_info(const typename C3t3::Edge& edge,
      bool& update_v0,
      bool& update_v1,
      const C3t3& c3t3,
      CellSelector cell_selector)
    {
      typedef typename C3t3::Vertex_handle Vertex_handle;

      Vertex_handle v0 = edge.first->vertex(edge.second);
      Vertex_handle v1 = edge.first->vertex(edge.third);

      int dim0 = c3t3.in_dimension(v0);
      int dim1 = c3t3.in_dimension(v1);

      std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
      std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

      update_v0 = false;
      update_v1 = false;

      bool is_v0_on_hull = is_on_hull(v0, c3t3);
      bool is_v1_on_hull = is_on_hull(v1, c3t3);

      //Same type imaginary or inside vertices
      if (dim0 == 3 && dim1 == 3)
      {
        if (is_v0_on_hull && is_v1_on_hull)//both endvertices are on hull
        {
          if (is_on_hull(edge, c3t3)) //edge also is on hull
          {
            update_v0 = true;
            update_v1 = true;
          }
        }
        else
        {
          if (!is_v0_on_hull) //v0 not on hull
            update_v0 = true;
          if (!is_v1_on_hull) //v1 not on hull
            update_v1 = true;
        }
        return;
      }
      //Feature edge case
      if (nb_si_v0 > 2 && nb_si_v1 > 2)
      {
        if (c3t3.is_in_complex(edge))
        {
          if (!topology_test(edge, c3t3, cell_selector))
            return;

          if (nb_si_v0 > nb_si_v1) {
            update_v1 = true;
          }
          else if (nb_si_v1 > nb_si_v0) {
            update_v0 = true;
          }
          else {
            update_v0 = true;
            update_v1 = true;
          }
        }
        return;
      }

      if (dim0 == 2 && dim1 == 2)
      {
        if (is_boundary(c3t3, edge, cell_selector))
        {
          if (!topology_test(edge, c3t3, cell_selector))
            return;
          Subdomain_relation subdomain_rel = compare_subdomains(v0, v1, c3t3);

          //Vertices on the same surface
          if (subdomain_rel == INCLUDES) {
            update_v1 = true;
          }
          else if (subdomain_rel == INCLUDED) {
            update_v0 = true;
          }
          else if (subdomain_rel == EQUAL)
          {
            if (c3t3.number_of_edges() == 0)
            {
              update_v0 = true;
              update_v1 = true;
            }
            else
            {
              bool v0_on_feature = is_on_feature(v0);
              bool v1_on_feature = is_on_feature(v1);

              if (v0_on_feature && v1_on_feature) {
                if (c3t3.is_in_complex(edge)) {
                  if (!c3t3.is_in_complex(v0))
                    update_v0 = true;
                  if (!c3t3.is_in_complex(v1))
                    update_v1 = true;
                }
              }
              else {
                if (!v0_on_feature) {
                  update_v0 = true;
                }
                if (!v1_on_feature) {
                  update_v1 = true;
                }
              }
            }
          }
        }

        return;
      }
      //In the case of mixte edges
      if (dim0 == 2 && dim1 == 3 && !is_v1_on_hull) {
        update_v1 = true;
        return;
      }

      if (dim1 == 2 && dim0 == 3 && !is_v0_on_hull) {
        update_v0 = true;
        return;
      }
    }


    template<typename C3T3>
    void print_subdomain_indices(const C3T3& c3t3)
    {
      typedef typename C3T3::Triangulation       Tr;
      typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;

      std::cout << "SUBDOMAINS : " << std::endl;
      unsigned int line_id = 0;
      for (Finite_cells_iterator cit = c3t3.triangulation().finite_cells_begin();
        cit != c3t3.triangulation().finite_cells_end();
        ++cit, ++line_id)
      {
        if (line_id % 10 == 0)
          std::cout << std::endl;
        std::cout << "\t" << cit->subdomain_index();
      }

    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    template <typename Bimap>
    void dump_edges(const Bimap& edges, const char* filename)
    {
      std::ofstream ofs(filename);
      ofs.precision(17);

      BOOST_FOREACH(typename Bimap::left_const_reference it, edges.left)
      {
        ofs << "2 " << it.first.first->point()
          << " " << it.first.second->point() << std::endl;
      }

      ofs.close();
    }
#endif
  }

  namespace internal
  {
    template<typename Tr, typename CellsSet, typename CellHandle>
    bool insert_in_cells(const CellHandle c, CellsSet& cells)
    {
      std::set<typename Tr::Vertex_handle> vertices;
      for (int i = 0; i < 4; ++i)
        vertices.insert(c->vertex(i));
      if (cells.find(vertices) != cells.end())
        return false;
      cells.insert(vertices);
      return true;
    }
  } // end internal

  namespace debug {
    // forward-declaration
    template<typename Tr, typename CellRange>
    void dump_cells(const CellRange& cells, const char* filename);
  }
  namespace helpers
  {

    template<typename K>
    void read_iso_cuboid(std::istream& is,
      CGAL::Iso_cuboid_3<K>& bbox)
    {
      typedef typename K::Point_3 Point_3;
      Point_3 p1, p2;
      double x, y, z;
      is >> x >> y >> z;
      p1 = Point_3(x, y, z);
      is >> x >> y >> z;
      p2 = Point_3(x, y, z);

      if (p1 < p2)
        bbox = CGAL::Iso_cuboid_3<K>(p1, p2);
      else
        bbox = CGAL::Iso_cuboid_3<K>(p2, p1);

      CGAL_assertion(p1 != p2);
    }

    template<typename Tr>
    void set_time_stamps(Tr& tr)
    {
      typedef typename Tr::Triangulation_data_structure::Vertex Vertex;
      typedef typename Tr::Triangulation_data_structure::Cell   Cell;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle   Cell_handle;

      CGAL::Time_stamper_impl<Vertex> v_ts;
      for (typename Tr::All_vertices_iterator vit = tr.all_vertices_begin();
        vit != tr.all_vertices_end();
        ++vit)
      {
        Vertex_handle vh = vit;
        Vertex* pv = &*vh;
        v_ts.initialize_time_stamp(pv);
        v_ts.set_time_stamp(pv);
      }
      CGAL::Time_stamper_impl<Cell> c_ts;
      for (typename Tr::All_cells_iterator cit = tr.all_cells_begin();
        cit != tr.all_cells_end();
        ++cit)
      {
        Cell_handle ch = cit;
        Cell* pc = &*ch;
        c_ts.initialize_time_stamp(pc);
        c_ts.set_time_stamp(pc);
      }
    }

    template<typename TDS_src, typename TDS_tgt>
    struct Vertex_converter
    {
      //This operator is used to create the vertex from v_src.
      typename TDS_tgt::Vertex operator()(const typename TDS_src::Vertex& v_src) const
      {
        typedef typename CGAL::Kernel_traits<
          typename TDS_src::Vertex::Point>::Kernel GT_src;
        typedef typename CGAL::Kernel_traits<
          typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
        CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

        typename TDS_tgt::Vertex v_tgt;
        v_tgt.set_point(conv(v_src.point()));
        v_tgt.set_time_stamp(-1);
        v_tgt.set_dimension(v_src.info());//-1 if unset, 0,1,2, or 3 if set
        return v_tgt;
      }
      //This operator is meant to be used in case heavy data should transferred to v_tgt.
      void operator()(const typename TDS_src::Vertex& v_src,
        typename TDS_tgt::Vertex& v_tgt) const
      {
        typedef typename CGAL::Kernel_traits<
          typename TDS_src::Vertex::Point>::Kernel GT_src;
        typedef typename CGAL::Kernel_traits<
          typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
        CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

        v_tgt.set_point(conv(v_src.point()));
        v_tgt.set_dimension(v_src.info());
      }
    };

    template<typename TDS_src, typename TDS_tgt>
    struct Cell_converter
    {
      //This operator is used to create the cell from c_src.
      typename TDS_tgt::Cell operator()(const typename TDS_src::Cell& c_src) const
      {
        typename TDS_tgt::Cell c_tgt;
        c_tgt.info() = c_src.info();
        c_tgt.input_cell() = c_src;
        c_tgt.set_time_stamp(-1);
        return c_tgt;
      }
      //This operator is meant to be used in case heavy data should transferred to c_tgt.
      void operator()(const typename TDS_src::Cell& c_src,
        typename TDS_tgt::Cell& c_tgt) const
      {
        c_tgt.info() = c_src.info();
        c_tgt.input_cell() = c_src;
      }
    };

    template<typename Tr, typename CellsSet, typename K>
    bool check_size_of_padding_box(const CellsSet& inside_cells,
      const CGAL::Iso_cuboid_3<K>& cuboid)
    {
      // all cells that are in intersecting_cells AND inside_cells
      // should NOT be clipped
      typedef typename CellsSet::value_type Cell_handle;

#ifdef CGAL_LIMITED_APERTURE_DEBUG
      std::vector<Cell_handle> cells;
#endif
      for (typename CellsSet::iterator cit = inside_cells.begin();
        cit != inside_cells.end(); ++cit)
      {
        Cell_handle c = *cit;
        for (int i = 0; i < 4; ++i)
        {
          if (cuboid.has_on_unbounded_side(c->vertex(i)->point()))
          {
#ifdef CGAL_LIMITED_APERTURE_DEBUG
            cells.push_back(c);
            break;
#else
            return false;
#endif
          }
        }
      }
#ifdef CGAL_LIMITED_APERTURE_DEBUG
      debug::dump_cells<Tr>(cells, "cells_from_padding_zone.mesh");
      return cells.empty();
#else
      return true;
#endif
    }

#ifdef CGAL_LIMITED_APERTURE_EDGE_SELECTION

    template<typename C3T3, typename CellCirculator>
    bool outer_box_criterion(const C3T3& c3t3,
      CellCirculator circ,
      CellCirculator end,
      const typename C3T3::Subdomain_index& imaginary_index)
    {
      std::size_t nb_imaginary = 0;
      std::size_t nb_total = 0;
      std::size_t nb_padding = 0;
      std::size_t nb_outside = 0;
      do
      {
        if (circ->subdomain_index() == imaginary_index)
          ++nb_imaginary;
        else if (!c3t3.is_in_complex(circ))
          ++nb_outside;
        else if (circ->info().padding())
          ++nb_padding;

        ++nb_total;
      } while (++circ != end);

      return nb_padding > 0 && (nb_imaginary + nb_outside) < nb_total;
    }

    template <typename C3T3>
    bool is_on_the_outer_box(const typename C3T3::Edge& e,
      const C3T3& c3t3,
      const typename C3T3::Subdomain_index& imaginary_index)
    {
      typedef typename C3T3::Triangulation::Cell_circulator Cell_circulator;
      Cell_circulator circ = c3t3.triangulation().incident_cells(e);
      Cell_circulator end = circ;

      return outer_box_criterion(c3t3, circ, end, imaginary_index);
    }

    template <typename C3T3>
    bool is_on_the_outer_box(const typename C3T3::Vertex_handle& v,
      const C3T3& c3t3,
      const typename C3T3::Subdomain_index& imaginary_index)
    {
      std::vector<typename C3T3::Cell_handle> cells;
      c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

      return outer_box_criterion(c3t3, cells.begin(), cells.end(), imaginary_index);
    }
#endif CGAL_LIMITED_APERTURE_EDGE_SELECTION

  }//end namespace helpers


  namespace debug
  {
    template<typename Tr>
    void rebuild_with_insert(Tr& tr, const int nbv_max)
    {
      typedef typename Tr::Point Point;
      std::vector<Point> points(tr.number_of_vertices());
      int index = 0;
      for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end(); ++vit)
      {
        points[index++] = vit->point();
      }

      tr.clear();
      for (int i = 0; i < nbv_max; ++i)
        tr.insert(points[i]);
    }

    template<typename Tr>
    void check_validity(const Tr& tr)
    {
      CGAL_assertion(tr.is_valid(true));
      for (typename Tr::All_vertices_iterator vit = tr.all_vertices_begin();
        vit != tr.all_vertices_end();
        ++vit)
      {
        typename Tr::Cell_handle c = vit->cell();
        CGAL_assertion(c->has_vertex(vit));
      }

      std::ofstream ofs("extra_cells.polylines.txt");
      std::set<typename Tr::Cell_handle> extra_cells;
      std::set<std::set<typename Tr::Vertex_handle> > cells;
      for (typename Tr::All_cells_iterator cit = tr.all_cells_begin();
        cit != tr.all_cells_end();
        ++cit)
      {
        typename Tr::Cell_handle c = cit;
        for (int i = 0; i < 4; ++i)
        {
          typename Tr::Cell_handle ci = c->neighbor(i);
          int j;
          CGAL_assertion(c->has_neighbor(ci, j));
          CGAL_assertion(i == j);
          j = ci->index(c);
          CGAL_assertion(ci->neighbor(j) == c);
          CGAL_assertion(ci->has_neighbor(c, j));
        }
        if (!internal::insert_in_cells<Tr>(c, cells))
        {
          extra_cells.insert(c);
          for (int j = 0; j < 4; ++j)
            dump_facet(std::make_pair(c, j), ofs);
        }
      }
      ofs.close();
    }

    template<typename ClippedCellsMap, typename Tr>
    void debug_infinite_facets(const ClippedCellsMap& clipped_cells,
      const Tr& tr)
    {
      typedef typename Tr::Point Point;
      //collect convex hull of tr edges (as pairs of ordered points)
      std::vector<std::pair<Point, Point>/*ordered pair*/> ch_edges;

      for (typename ClippedCellsMap::const_iterator cmit = clipped_cells.begin();
        cmit != clipped_cells.end();
        ++cmit)
      {
        typedef typename ClippedCellsMap::mapped_type CellTr;
        const CellTr& ctr = cmit->second;

        typename ClippedCellsMap::key_type cell = cmit->first;

        for (typename CellTr::Finite_edges_iterator eit = ctr.finite_edges_begin();
          eit != ctr.finite_edges_end();
          ++eit)
        {
          Point p1 = (eit->first)->vertex(eit->second)->point();
          Point p2 = (eit->first)->vertex(eit->third)->point();
          if (p2 < p1)
            std::swap(p1, p2); //make sure that p1 <= p2

          int vi = 0, vj = 0;
          for (; vi < 4; ++vi)
          {
            if (cell->vertex(vi)->point() == p1)
              break;
          }
          if (vi == 4)
            continue;
          for (; vj < 4; ++vj)
          {
            if (cell->vertex(vj)->point() == p2)
              break;
          }
          if (vj == 4)
            continue;

          int vk = Tr::next_around_edge(vi, vj);
          int vl = Tr::next_around_edge(vj, vi);
          if (tr.is_infinite(cell->neighbor(vk)))
            ch_edges.push_back(std::make_pair(p1, p2));
          if (tr.is_infinite(cell->neighbor(vl)))
            ch_edges.push_back(std::make_pair(p1, p2));
        }
      }

      //check that each edge appears exactly twice
      std::sort(ch_edges.begin(), ch_edges.end());
      bool twice_each = (ch_edges.size() % 2 == 0);
      for (std::size_t i = 0; i < ch_edges.size() - 1; i = i + 2)
      {
        if (ch_edges[i] != ch_edges[i + 1])
          twice_each = false;
      }
      if (!twice_each)
      {
        for (std::size_t i = 0; i < ch_edges.size(); ++i)
        {
          std::cout << i << "\t"
            << ch_edges[i].first << " " << ch_edges[i].second
            << std::endl;
        }
      }
      CGAL_assertion(twice_each);
    }

    template<typename Tr>
    void dump_surface_off(const Tr& tr, const char* filename)
    {
      typedef typename Tr::Vertex_handle              Vertex_handle;
      typedef typename Tr::Cell_handle                Cell_handle;
      typedef typename Tr::Finite_facets_iterator     Finite_facets_iterator;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      //collect vertices
      Bimap_t vertices;
      std::size_t nbf = 0;
      int index = 0;
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        if (tr.is_infinite(c) || tr.is_infinite(c->neighbor(i)))
        {
          nbf++;
          for (int j = 1; j < 4; ++j)
          {
            Vertex_handle vij = c->vertex((i + j) % 4);
            if (vertices.left.find(vij) == vertices.left.end())
              vertices.left.insert(value_type(vij, index++));
          }
        }
      }

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << vertices.left.size() << " " << nbf << " 0" << std::endl << std::endl;

      // write vertices
      for (typename Bimap_t::right_iterator vit = vertices.right.begin();
        vit != vertices.right.end(); ++vit)
      {
        ofs << vit->second->point() << std::endl;
      }

      //write facets
      std::size_t nbf_print = 0;
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        if (tr.is_infinite(c) || tr.is_infinite(c->neighbor(i)))
        {
          ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
            << vertices.left.at(c->vertex((i + 2) % 4)) << " "
            << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
          ++nbf_print;
        }
      }
      CGAL_assertion(nbf == nbf_print);

      ofs.close();
    }

    template<typename Tr>
    void dump_cells_off(const Tr& tr, const char* filename)
    {
      typedef typename Tr::Vertex_handle              Vertex_handle;
      typedef typename Tr::Cell_handle                Cell_handle;
      typedef typename Tr::Finite_facets_iterator     Finite_facets_iterator;
      typedef typename Tr::Finite_vertices_iterator   Finite_vertices_iterator;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << tr.number_of_vertices()
        << " " << tr.number_of_finite_facets() << " 0" << std::endl << std::endl;

      //collect and write vertices
      Bimap_t vertices;
      int index = 0;
      for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end(); ++vit)
      {
        vertices.left.insert(value_type(vit, index++));
        ofs << vit->point().x() << " "
          << vit->point().y() << " "
          << vit->point().z() << std::endl;
      }

      //write facets
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
          << vertices.left.at(c->vertex((i + 2) % 4)) << " "
          << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
      }
      ofs.close();
    }

    template<typename Tr, typename CellRange, typename IndexRange>
    void dump_cells(const CellRange& cells,
      const IndexRange& indices,
      const char* filename)
    {
      typedef typename Tr::Vertex_handle                         Vertex_handle;
      typedef typename Tr::Point                                 Point;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      CGAL_assertion(indices.empty() || cells.size() == indices.size());

      //collect vertices
      Bimap_t vertices;
      int index = 1;
      for (typename CellRange::const_iterator cit = cells.begin();
        cit != cells.end();
        ++cit)
      {
        for (int i = 0; i < 4; ++i)
        {
          Vertex_handle vi = (*cit)->vertex(i);
          if (vertices.left.find(vi) == vertices.left.end())
            vertices.left.insert(value_type(vi, index++));
        }
      }

      //write cells
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "MeshVersionFormatted 1" << std::endl;
      ofs << "Dimension 3" << std::endl;
      ofs << "Vertices" << std::endl << vertices.size() << std::endl;
      for (typename Bimap_t::right_const_iterator vit = vertices.right.begin();
        vit != vertices.right.end();
        ++vit)
      {
        const Point& p = vit->second->point();
        ofs << p.x() << " " << p.y() << " " << p.z() << " 2" << std::endl;
      }
      ofs << "Tetrahedra " << std::endl << cells.size() << std::endl;
      typename IndexRange::const_iterator iit = indices.begin();
      for (typename CellRange::const_iterator cit = cells.begin();
        cit != cells.end();
        ++cit)
      {
        ofs << vertices.left.at((*cit)->vertex(0))
          << " " << vertices.left.at((*cit)->vertex(1))
          << " " << vertices.left.at((*cit)->vertex(2))
          << " " << vertices.left.at((*cit)->vertex(3));

        if (iit == indices.end())
          ofs << " 1" << std::endl;
        else
        {
          // std::cerr << "Cell #" << (cit - cells.begin())
          //           << " has original index " << *iit << std::endl;
          ofs << " " << (*iit) << std::endl;
          ++iit;
        }
      }
      ofs << "End" << std::endl;
      ofs.close();
    }

    template<typename Tr, typename CellRange>
    void dump_cells(const CellRange& cells, const char* filename)
    {
      std::vector<int> indices;
      dump_cells<Tr>(cells, indices, filename);
    }

    template<typename Tr>
    void dump_cells_in_complex(const Tr& tr, const char* filename)
    {
      std::vector<typename Tr::Cell_handle> cells;
      std::vector<typename Tr::Cell::Subdomain_index> indices;

      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->subdomain_index() > 0)
        {
          cells.push_back(cit);
          indices.push_back(cit->subdomain_index());
        }
      }
      dump_cells<Tr>(cells, indices, filename);
    }

    template<typename C3t3>
    void dump_facets_in_complex(const C3t3& c3t3, const char* filename)
    {
      typedef typename C3t3::Triangulation              Tr;
      typedef typename Tr::Vertex_handle                Vertex_handle;
      typedef typename Tr::Cell_handle                  Cell_handle;
      typedef typename C3t3::Facets_in_complex_iterator Facets_in_complex_iterator;
      typedef boost::bimap<Vertex_handle, int>          Bimap_t;
      typedef typename Bimap_t::left_map::value_type  value_type;

      //collect vertices
      Bimap_t vertices;
      std::size_t nbf = 0;
      int index = 0;
      for (Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin();
        fit != c3t3.facets_in_complex_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;

        nbf++;
        for (int j = 1; j < 4; ++j)
        {
          Vertex_handle vij = c->vertex((i + j) % 4);
          if (vertices.left.find(vij) == vertices.left.end())
            vertices.left.insert(value_type(vij, index++));
        }
      }

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << vertices.left.size() << " " << nbf << " 0" << std::endl << std::endl;

      // write vertices
      for (typename Bimap_t::right_iterator vit = vertices.right.begin();
        vit != vertices.right.end(); ++vit)
      {
        ofs << vit->second->point() << std::endl;
      }

      //write facets
      std::size_t nbf_print = 0;
      for (Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin();
        fit != c3t3.facets_in_complex_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
          << vertices.left.at(c->vertex((i + 2) % 4)) << " "
          << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
        ++nbf_print;
      }
      CGAL_assertion(nbf == nbf_print);

      ofs.close();
    }

    template<typename C3T3>
    void dump_edges_in_complex(const C3T3& c3t3, const char* filename)
    {
      std::ofstream ofs(filename);
      ofs.precision(17);
      for (typename C3T3::Edges_in_complex_iterator eit = c3t3.edges_in_complex_begin();
        eit != c3t3.edges_in_complex_end(); ++eit)
      {
        const typename C3T3::Edge& e = *eit;
        ofs << "2 "
          << e.first->vertex(e.second)->point() << " "
          << e.first->vertex(e.third)->point() << "\n";
      }
      ofs.close();
    }

    template<typename Tr>
    void dump_vertices_by_dimension(const Tr& tr, const char* prefix)
    {
      typedef typename Tr::Vertex_handle Vertex_handle;
      std::vector< std::vector<Vertex_handle> > vertices_per_dimension(4);

      for (typename Tr::Finite_vertices_iterator
        vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end();
        ++vit)
      {
        //vertices_per_dimension[vit->info()].push_back(vit);
        vertices_per_dimension[vit->in_dimension()].push_back(vit);
      }

      for (int i = 0; i < 4; ++i)
      {
        //dimension is i
        const std::vector<Vertex_handle>& vertices_di = vertices_per_dimension[i];

        std::cout << "Dimension " << i << " : " << vertices_di.size() << std::endl;

        std::ostringstream oss;
        oss << prefix << "_dimension_" << i << ".off";

        std::ofstream ofs(oss.str());
        ofs.precision(17);
        ofs << "OFF" << std::endl;
        ofs << vertices_di.size() << " 0 0" << std::endl << std::endl;

        for (std::size_t j = 0; j < vertices_di.size(); ++j)
        {
          ofs << vertices_di[j]->point() << std::endl;
        }

        ofs.close();
      }
    }

    template<typename Tr>
    void dump_triangulation_cells(const Tr& tr, const char* filename)
    {
      std::vector<typename Tr::Cell_handle> cells(tr.number_of_finite_cells());
      int i = 0;
      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        cells[i++] = cit;
      }
      dump_cells<Tr>(cells, filename);
    }

    template<typename Tr>
    void dump_without_imaginary(const Tr& tr, const char* filename,
      const int imaginary_index)
    {
      std::vector<typename Tr::Cell_handle> cells;
      std::vector<std::ptrdiff_t> indices;

      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->subdomain_index() > 0
          && cit->subdomain_index() != imaginary_index)
        {
          cells.push_back(cit);
          indices.push_back(1);
            //cit->info().padding() ?
            //-1 :
            //cit->info().original_index());
        }
      }
      dump_cells<Tr>(cells, indices, filename);
    }

    template<typename Tr>
    void dump_padding_cells(const Tr& tr, const char* filename)
    {
      std::vector<typename Tr::Cell_handle> cells;
      std::vector<typename Tr::Cell::Subdomain_index> indices;

      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->info().padding())
        {
          cells.push_back(cit);
          if (cit->subdomain_index() > 0)
            indices.push_back(cit->subdomain_index());
          else
            indices.push_back(1);
        }
      }
      dump_cells<Tr>(cells, indices, filename);
    }

    template <typename Tr, typename Isocuboid>
    void dump_non_padding_plus_the_outer_bbox(const Tr& tr,
      const Isocuboid& bbox,
      const char* filename)
    {
      typedef typename Tr::Vertex_handle                         Vertex_handle;
      typedef typename Tr::Point                                 Point;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      Bimap_t vertices;
      int index = 9; // because we output first the 8 vertices of the
                     // bbox
      std::size_t nb_of_cells = 0;
      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->info().padding()) {
          continue;
        }
        ++nb_of_cells;
        for (int i = 0; i < 4; ++i)
        {
          Vertex_handle vi = cit->vertex(i);
          if (vertices.left.find(vi) == vertices.left.end())
            vertices.left.insert(value_type(vi, index++));
        }
      }
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "MeshVersionFormatted 1\n"
        << "Dimension 3\n"
        << "Vertices\n"
        << vertices.size() + 8 << std::endl;
      const CGAL::cpp11::array<int, 8> indices = { 0, 3, 2, 1, 5, 4, 7, 6 };
      for (int i = 0; i < 8; ++i) {
        const typename Tr::Point_3 p = bbox[indices[i]];
        ofs << p.x() << " " << p.y() << " " << p.z() << " 1" << std::endl;
      }
      for (typename Bimap_t::right_const_iterator vit = vertices.right.begin();
        vit != vertices.right.end();
        ++vit)
      {
        const Point& p = vit->second->point();
        ofs << p.x() << " " << p.y() << " " << p.z() << " 2" << std::endl;
      }
      ofs << "Triangles\n"
        << "12\n"
        << "1 2 4 1\n"
        << "4 2 3 1\n"
        << "1 5 2 1\n"
        << "2 5 6 1\n"
        << "4 3 8 1\n"
        << "8 3 7 1\n"
        << "5 1 4 1\n"
        << "8 5 4 1\n"
        << "7 5 8 1\n"
        << "7 6 5 1\n"
        << "2 6 7 1\n"
        << "3 2 7 1\n";
      ofs << "Tetrahedra" << std::endl << nb_of_cells << std::endl;
      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->info().padding()) continue;
        ofs << vertices.left.at(cit->vertex(0)) << " "
          << vertices.left.at(cit->vertex(1)) << " "
          << vertices.left.at(cit->vertex(2)) << " "
          << vertices.left.at(cit->vertex(3)) << " "
          << cit->info().original_index() << std::endl;
      }
      ofs << "End" << std::endl;
      ofs.close();
    }


    template<typename VertexPairsSet>
    void dump_edges(const VertexPairsSet& edges, const char* filename)
    {
      std::ofstream ofs(filename);
      BOOST_FOREACH(typename VertexPairsSet::key_type vp, edges)
      {
        ofs << "2 " << vp.first->point()
          << " " << vp.second->point() << std::endl;
      }
      ofs.close();
    }
  }// end namespace debug

}
}

#endif //CGAL_INTERNAL_TET_REMESHING_HELPERS_H