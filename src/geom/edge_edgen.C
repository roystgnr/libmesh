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
#include "libmesh/edge_edgen.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

// ------------------------------------------------------------
// EdgeN class static member initializations
template <unsigned int N>
const int EdgeN::num_nodes;


template <unsigned int N>
bool EdgeN::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}


template <unsigned int N>
bool EdgeN::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}


template <unsigned int N>
bool EdgeN::is_face(const unsigned int ) const
{
  return false;
}


template <unsigned int N>
bool EdgeN::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, num_nodes);
  return (s == n);
}


template <unsigned int N>
bool EdgeN::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}


template <unsigned int N>
bool EdgeN::has_affine_map() const
{
  const Point vxi = (this->point(1) - this->point(0))/(N-1);
  for (auto i : make_range(1,N-1))
    {
      const Point v = vxi * i;
      if (!v.relative_fuzzy_equals
            (this->point(i+1) - this->point(0), affine_tol))
        return false;
    }

  return true;
}


template <unsigned int N>
Order EdgeN::default_order() const
{
  return static_cast<Order>(N-1);
}


template <unsigned int N>
void EdgeN::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_not_implemented();
}


template <unsigned int N>
std::pair<unsigned short int, unsigned short int>
EdgeN::second_order_child_vertex (const unsigned int) const
{
  libmesh_not_implemented();
  return std::pair<unsigned short int, unsigned short int>(0,0);
}


template <unsigned int N>
BoundingBox EdgeN::loose_bounding_box () const
{
  libmesh_not_implemented();
  return BoundingBox(0, 0);
}


template <unsigned int N>
void EdgeN::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  for (auto j : make_range((N-2)/2))
    swap2nodes(2+i, N-i);
  swap2neighbors(0,1);
  swap2boundarysides(0,1,boundary_info);
}


} // namespace libMesh
