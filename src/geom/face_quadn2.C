// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_quad9.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// QuadNSq class static member initializations
template <unsigned int N>
const int QuadNSq::num_nodes;

template <unsigned int N>
const int QuadNSq::nodes_per_side;


// ------------------------------------------------------------
// QuadNSq class member functions

template <unsigned int N>
bool QuadNSq::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

template <unsigned int N>
bool QuadNSq::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  if (i >= 4*N-4)
    return false;
  return true;
}

template <unsigned int N>
bool QuadNSq::is_face(const unsigned int i) const
{
  if (i >= 4*N-4)
    return true;
  return false;
}

template <unsigned int N>
bool QuadNSq::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  if ((n-4)/(N-2) == s)
    return true;
  return false;
}


template <unsigned int N>
std::vector<unsigned int>
QuadNSq::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  std::vector<unsigned int> returnval(N);
  returnval[0] = s;
  returnval[1] = (s+1)%4;
  std::iota(returnval.begin()+2, returnval.end(), 4+s*(N-2));
  return returnval;
}


template <unsigned int N>
std::vector<unsigned int>
QuadNSq::nodes_on_edge(const unsigned int e) const
{
  return this->nodes_on_side(e);
}


template <unsigned int N>
bool QuadNSq::has_affine_map() const
{
  // make sure corners form a parallelogram
  const Point v01 = this->point(1) - this->point(0);
  const Point v03 = this->point(3) - this->point(0);
  if (!v01.relative_fuzzy_equals(this->point(2) - this->point(3), affine_tol))
    return false;

  const Point vxi = v01 / (N-1);
  const Point veta = v03 / (N-1);

  // make sure sides are straight
  for (auto i : make_range(1,N-1))
    {
      const Point vhoriz = vxi * i;
      const Point vvert = veta * i;
      if (!vhoriz.relative_fuzzy_equals(this->point(3+i) - this->point(0), affine_tol))
        return false;
      if (!vhoriz.relative_fuzzy_equals(this->point(3*N-2-i) - this->point(3), affine_tol))
        return false;
      if (!vvert.relative_fuzzy_equals(this->point(N+1+i) - this->point(1), affine_tol))
        return false;
      if (!vvert.relative_fuzzy_equals(this->point(4*N-4-i) - this->point(0), affine_tol))
        return false;

      // Make sure the interior isn't distorted
      for (auto j : make_range(1,N-1))
        {
          const Point v = vvert + vxi*j;
          if (!v.relative_fuzzy_equals(this->point(3*N-3+(N-2)*i+j) - this->point(0), affine_tol))
            return false;
        }
    }

  return true;
}



template <unsigned int N>
Order QuadNSq::default_order() const
{
  return static_cast<Order>(N-1);
}



template <unsigned int N>
unsigned int QuadNSq::local_side_node(unsigned int side,
                                      unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, nodes_per_side);

  if (side_node == 0)
    return side;

  if (side_node == 1)
    return (side+1)%4;

  return side * (N-2) + 2;
}



template <unsigned int N>
std::unique_ptr<Elem> QuadNSq::build_side_ptr (const unsigned int i,
                                               bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  libmesh_assert(!proxy); // deprecated long ago

  std::unique_ptr<Elem> side = Elem::build(this->side_type(i));
  this->build_side_ptr(side, i);
  return side;
}



template <unsigned int N>
void QuadNSq::build_side_ptr (std::unique_ptr<Elem> & side,
                              const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  const ElemType edgeNtype = this->side_type(i);

  if (side->type() != edgeNtype)
    {
      side = this->build_side_ptr(i);
      return;
    }

  for (auto n : side->node_index_range())
    side->set_node(n) = this->node_ptr(this->local_side_node(i,n));

  side->set_parent(nullptr);
  side->set_interior_parent(this);

  side->set_mapping_type(this->mapping_type());
  side->subdomain_id() = this->subdomain_id();
#ifdef LIBMESH_ENABLE_AMR
  side->set_p_level(this->p_level());
#endif
}



template <unsigned int N>
void QuadNSq::connectivity(const unsigned int sf,
                           const IOPackage iop,
                           std::vector<dof_id_type> & conn) const
{
  libmesh_not_implemented();
}



template <unsigned int N>
BoundingBox QuadNSq::loose_bounding_box () const
{
  // This might have curved edges, or might be a curved surface in
  // 3-space, in which case the full bounding box can be larger than
  // the bounding box of just the nodes.
  //
  //
  // FIXME
  libmesh_not_implemented();
  return BoundingBox(0, 0);
}



template <unsigned int N>
unsigned int QuadNSq::n_second_order_adjacent_vertices (const unsigned int n) const
{
  libmesh_not_implemented();
}



template <unsigned int N>
unsigned short int QuadNSq::second_order_adjacent_vertex (const unsigned int n,
                                                          const unsigned int v) const
{
  libmesh_not_implemented();
}



template <unsigned int N>
std::pair<unsigned short int, unsigned short int>
QuadNSq::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_not_implemented();
}


template <unsigned int N>
void QuadNSq::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 4);

  for (auto i : make_range(perm_num))
    {
      swap4nodes(0,1,2,3);
      for (auto j : make_range(N-2))
        {
          swap4nodes(4+j,4+j+N-2,4+j+2*N-4,4+j+3*N-6);
        }
      swap4neighbors(0,1,2,3);
    }
}


template <unsigned int N>
void QuadNSq::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  for (auto i : make_range(N-2))
    swap2nodes(N+1+i, 4*(N-1)-i);

  for (auto j : make_range((N-2)/2))
    {
      swap2nodes(3+i, N+2-i);
      swap2nodes(2*N-1+i, 3*N-2-i);
    }

  swap2neighbors(1,3);
  swap2boundarysides(1,3,boundary_info);
  swap2boundaryedges(1,3,boundary_info);
}


template <unsigned int N>
unsigned int QuadNSq::center_node_on_side(const unsigned short side) const
{
  libmesh_not_implemented();
  return invalid_uint;
}


template <unsigned int N>
ElemType QuadNSq::side_type (const unsigned int libmesh_dbg_var(s)) const
{
  libmesh_assert_less (s, 4);
  if (N == 2)
    return EDGE2;
  if (N == 3)
    return EDGE3;
  if (N == 4)
    return EDGE4;

  return static_cast<ElemType>(static_cast<unsigned int>(EDGEN)+N);
}


} // namespace libMesh
