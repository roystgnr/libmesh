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



#ifndef LIBMESH_EDGE_EDGEN_H
#define LIBMESH_EDGE_EDGEN_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/edge.h"

namespace libMesh
{

/**
 * The \p EdgeN is an element in 1D composed of N nodes. It is numbered
 * e.g. like this for 5 nodes:
 *
 * \verbatim
 *   EDGEN: o--o--o--o--o        o---> xi
 *          0  2  3  4  1
 * \endverbatim
 *
 * xi in [-1,1] is the reference element coordinate associated with
 * the given numbering.
 *
 * \author Roy H. Stogner
 * \date 2025
 * \brief A 1D geometric element with N nodes.
 */
template <unsigned int N>
class EdgeN : public Edge
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Edge3 (Elem * p=nullptr) :
    Edge(num_nodes, p, _nodelinks_data) {}

  Edge3 (Edge3 &&) = delete;
  Edge3 (const Edge3 &) = delete;
  Edge3 & operator= (const Edge3 &) = delete;
  Edge3 & operator= (Edge3 &&) = delete;
  virtual ~Edge3() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override
  {
    libmesh_assert_less(i, this->n_nodes());
    if (i < 2)
      return Point(2.0f*Real(i)-1.0f,0,0);
    return Point(0,0,0);
  }

  /**
   * \returns 3.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 2.
   */
  virtual unsigned int n_sub_elem() const override { return 2; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (always true in 1D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns the appropriate enum depending on N
   */
  virtual ElemType type() const override;

  /**
   * \returns N-1
   */
  virtual Order default_order() const override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { libmesh_not_implemented(); }

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int,
                                                           const unsigned int v) const override
  { libmesh_not_implemented(); return 0; }

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the edge.
   */
  virtual BoundingBox loose_bounding_box () const override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p false.  This is a finite element.
   */
  virtual bool infinite () const override { return false; }

#endif

  static const int num_nodes = N;

  virtual void flip(BoundaryInfo *) override final;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Matrix used to create the elements children.
   */
  virtual Real embedding_matrix (const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int k) const override
  { libmesh_not_implemented(); return 0; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

} // namespace libMesh


#endif // LIBMESH_EDGE_EDGEN_H
