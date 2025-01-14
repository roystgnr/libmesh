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



#ifndef LIBMESH_FACE_QUADNSQ_H
#define LIBMESH_FACE_QUADNSQ_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_quad.h"

namespace libMesh
{

/**
 * The \p QUADNSQ<N> is a quadrilateral element in 2D composed of N^2 nodes.
 * It is numbered e.g. like this for N=4
 * \verbatim
 *           3   9   8   2
 *  QUADNSQ: o---o---o---o
 *           |           |           eta
 *        10 o 14o 15o   o 7          ^
 *           |           |            |
 *        11 o 12o 13o   o 6          |
 *           |           |            o---> xi
 *           o---o---o---o
 *           0   4   5   1
 * \endverbatim
 * (xi, eta) in [-1,1]^2 are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Roy H Stogner
 * \date 2025
 * \brief A 2D quadrilateral element with 16 nodes.
 */
class QuadNSq : public Quad
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  QuadNSq (Elem * p=nullptr) :
    Quad(num_nodes, p, _nodelinks_data) {}

  QuadNSq (QuadNSq &&) = delete;
  QuadNSq (const QuadNSq &) = delete;
  QuadNSq & operator= (const QuadNSq &) = delete;
  QuadNSq & operator= (QuadNSq &&) = delete;
  virtual ~QuadNSq() = default;

  /**
   * \returns the appropriate enum depending on N
   */
  virtual ElemType type () const override;

  /**
   * \returns N^2.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_sub_elem() const override { return 4; }

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

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns N-1
   */
  virtual Order default_order() const override;

  /**
   * \returns The local node id for node \p side_node on side \p side of
   * this Elem.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=false) override;

  /**
   * Rebuilds an EDGEN coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for edge nodes and 4 for the face node.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const override;

  /**
   * \returns The element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 4 \le n < N^2 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element. See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * Geometric constants for QuadNSq.
   */
  static const int num_nodes = N*N;
  static const int nodes_per_side = N;

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const override;

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

  unsigned int center_node_on_side(const unsigned short side) const override final;

  ElemType side_type (const unsigned int s) const override final;

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


#endif // LIBMESH_FACE_QUADNSQ_H
