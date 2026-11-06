// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe.h"
#include "libmesh/inf_fe.h"
#include "libmesh/libmesh_logging.h"

// For projection code:
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/int_range.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/tensor_value.h"
#include "libmesh/threads.h"
#include "libmesh/fe_type.h"
#include "libmesh/enum_to_string.h"

// C++ Includes
#include <memory>

namespace libMesh
{



// ------------------------------------------------------------
// FEBase class members
template <>
std::unique_ptr<FEGenericBase<Real>>
FEGenericBase<Real>::build (const unsigned int dim,
                            const FEType & fet)
{
  switch (dim)
    {
      // 0D
    case 0:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<0,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<0,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<0,LAGRANGE>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<0,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_GLL:
            return std::make_unique<FE<0,L2_LAGRANGE_GLL>>(fet);

          case LAGRANGE_GLL:
            return std::make_unique<FE<0,LAGRANGE_GLL>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<0,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<0,L2_HIERARCHIC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<0,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<0,MONOMIAL>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<0,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<0,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<0,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<0>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<0>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }
      // 1D
    case 1:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<1,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<1,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<1,LAGRANGE>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<1,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_GLL:
            return std::make_unique<FE<1,L2_LAGRANGE_GLL>>(fet);

          case LAGRANGE_GLL:
            return std::make_unique<FE<1,LAGRANGE_GLL>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<1,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<1,L2_HIERARCHIC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<1,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<1,MONOMIAL>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<1,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<1,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<1,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<1>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<1>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }


      // 2D
    case 2:
      {
        switch (fet.family)
          {
          case CLOUGH:
            return std::make_unique<FE<2,CLOUGH>>(fet);

          case HERMITE:
            return std::make_unique<FE<2,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<2,LAGRANGE>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<2,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_GLL:
            return std::make_unique<FE<2,L2_LAGRANGE_GLL>>(fet);

          case LAGRANGE_GLL:
            return std::make_unique<FE<2,LAGRANGE_GLL>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<2,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<2,L2_HIERARCHIC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<2,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<2,MONOMIAL>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<2,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<2,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<2,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<2>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<2>>(fet);

          case SUBDIVISION:
            return std::make_unique<FESubdivision>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }


      // 3D
    case 3:
      {
        switch (fet.family)
          {
          case CLOUGH:
            libmesh_error_msg("ERROR: Clough-Tocher elements currently only support 1D and 2D");

          case HERMITE:
            return std::make_unique<FE<3,HERMITE>>(fet);

          case LAGRANGE:
            return std::make_unique<FE<3,LAGRANGE>>(fet);

          case L2_LAGRANGE:
            return std::make_unique<FE<3,L2_LAGRANGE>>(fet);

          case L2_LAGRANGE_GLL:
            return std::make_unique<FE<3,L2_LAGRANGE_GLL>>(fet);

          case LAGRANGE_GLL:
            return std::make_unique<FE<3,LAGRANGE_GLL>>(fet);

          case HIERARCHIC:
            return std::make_unique<FE<3,HIERARCHIC>>(fet);

          case L2_HIERARCHIC:
            return std::make_unique<FE<3,L2_HIERARCHIC>>(fet);

          case SIDE_HIERARCHIC:
            return std::make_unique<FE<3,SIDE_HIERARCHIC>>(fet);

          case MONOMIAL:
            return std::make_unique<FE<3,MONOMIAL>>(fet);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
          case SZABAB:
            return std::make_unique<FE<3,SZABAB>>(fet);

          case BERNSTEIN:
            return std::make_unique<FE<3,BERNSTEIN>>(fet);

          case RATIONAL_BERNSTEIN:
            return std::make_unique<FE<3,RATIONAL_BERNSTEIN>>(fet);
#endif

          case XYZ:
            return std::make_unique<FEXYZ<3>>(fet);

          case SCALAR:
            return std::make_unique<FEScalar<3>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }
}



template <>
std::unique_ptr<FEGenericBase<RealGradient>>
FEGenericBase<RealGradient>::build (const unsigned int dim,
                                    const FEType & fet)
{
  switch (dim)
    {
      // 0D
    case 0:
      {
        switch (fet.family)
          {
          case HIERARCHIC_VEC:
            return std::make_unique<FEHierarchicVec<0>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FEL2HierarchicVec<0>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FELagrangeVec<0>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FEL2LagrangeVec<0>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FEMonomialVec<0>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }
    case 1:
      {
        switch (fet.family)
          {
          case HIERARCHIC_VEC:
            return std::make_unique<FEHierarchicVec<1>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FEL2HierarchicVec<1>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FELagrangeVec<1>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FEL2LagrangeVec<1>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FEMonomialVec<1>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }
    case 2:
      {
        switch (fet.family)
          {
          case HIERARCHIC_VEC:
            return std::make_unique<FEHierarchicVec<2>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FEL2HierarchicVec<2>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FELagrangeVec<2>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FEL2LagrangeVec<2>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FEMonomialVec<2>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<2>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<2>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<2>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }
    case 3:
      {
        switch (fet.family)
          {
          case HIERARCHIC_VEC:
            return std::make_unique<FEHierarchicVec<3>>(fet);

          case L2_HIERARCHIC_VEC:
            return std::make_unique<FEL2HierarchicVec<3>>(fet);

          case LAGRANGE_VEC:
            return std::make_unique<FELagrangeVec<3>>(fet);

          case L2_LAGRANGE_VEC:
            return std::make_unique<FEL2LagrangeVec<3>>(fet);

          case MONOMIAL_VEC:
            return std::make_unique<FEMonomialVec<3>>(fet);

          case NEDELEC_ONE:
            return std::make_unique<FENedelecOne<3>>(fet);

          case RAVIART_THOMAS:
            return std::make_unique<FERaviartThomas<3>>(fet);

          case L2_RAVIART_THOMAS:
            return std::make_unique<FEL2RaviartThomas<3>>(fet);

          default:
            libmesh_error_msg("ERROR: Bad FEType.family == " << Utility::enum_to_string(fet.family));
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    } // switch(dim)
}







#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


template <>
std::unique_ptr<FEGenericBase<Real>>
FEGenericBase<Real>::build_InfFE (const unsigned int dim,
                                  const FEType & fet)
{
  switch (dim)
    {

      // 1D
    case 1:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Can't build an infinite element with FEFamily = " << Utility::enum_to_string(fet.radial_family));

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<1,JACOBI_20_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<1,JACOBI_30_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<1,LEGENDRE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<1,LAGRANGE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Can't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << Utility::enum_to_string(fet.radial_family));
          }
      }




      // 2D
    case 2:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Can't build an infinite element with FEFamily = " << Utility::enum_to_string(fet.radial_family));

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<2,JACOBI_20_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<2,JACOBI_30_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<2,LEGENDRE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<2,LAGRANGE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << Utility::enum_to_string(fet.radial_family));
          }
      }




      // 3D
    case 3:
      {
        switch (fet.radial_family)
          {
          case INFINITE_MAP:
            libmesh_error_msg("ERROR: Don't build an infinite element with FEFamily = " << Utility::enum_to_string(fet.radial_family));

          case JACOBI_20_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<3,JACOBI_20_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case JACOBI_30_00:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<3,JACOBI_30_00,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LEGENDRE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<3,LEGENDRE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          case LAGRANGE:
            {
              switch (fet.inf_map)
                {
                case CARTESIAN:
                  return std::make_unique<InfFE<3,LAGRANGE,CARTESIAN>>(fet);

                default:
                  libmesh_error_msg("ERROR: Don't build an infinite element with InfMapType = " << Utility::enum_to_string(fet.inf_map));
                }
            }

          default:
            libmesh_error_msg("ERROR: Bad FEType.radial_family= " << Utility::enum_to_string(fet.radial_family));
          }
      }

    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }
}



template <>
std::unique_ptr<FEGenericBase<RealGradient>>
FEGenericBase<RealGradient>::build_InfFE (const unsigned int,
                                          const FEType & )
{
  // No vector types defined... YET.
  libmesh_not_implemented();
  return std::unique_ptr<FEVectorBase>();
}

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


template <typename OutputType>
void FEGenericBase<OutputType>::compute_shape_functions (const Elem * elem,
                                                         const std::vector<Point> & qp)
{
  //-------------------------------------------------------------------------
  // Compute the shape function values (and derivatives)
  // at the Quadrature points.  Note that the actual values
  // have already been computed via init_shape_functions

  // Start logging the shape function computation
  LOG_ONE_ENTITY_SCOPE("compute_shape_functions()", "FE");

  this->determine_calculations();

  if (calculate_phi)
    this->_fe_trans->map_phi(this->dim, elem, qp, (*this), this->phi, this->_add_p_level_in_reinit);

  if (calculate_dphi)
    this->_fe_trans->map_dphi(this->dim, elem, qp, (*this), this->dphi,
                              this->dphidx, this->dphidy, this->dphidz);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (calculate_d2phi)
    this->_fe_trans->map_d2phi(this->dim, qp, (*this), this->d2phi,
                               this->d2phidx2, this->d2phidxdy, this->d2phidxdz,
                               this->d2phidy2, this->d2phidydz, this->d2phidz2);
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

  // Only compute curl for vector-valued elements
  if (calculate_curl_phi && TypesEqual<OutputType,RealGradient>::value)
    this->_fe_trans->map_curl(this->dim, elem, qp, (*this), this->curl_phi);

  // Only compute div for vector-valued elements
  if (calculate_div_phi && TypesEqual<OutputType,RealGradient>::value)
    this->_fe_trans->map_div(this->dim, elem, qp, (*this), this->div_phi);
}


// Here, we rely on the input \p phi_vals for accurate integration of the mass matrix.
// This is because in contact, we often have customized qrule for mortar segments
// due to deformation of the element, and the size of \p phi_vals for the secondary
// element changes accordingly.
template <>
void FEGenericBase<Real>::compute_dual_shape_coeffs (const std::vector<Real> & JxW, const std::vector<std::vector<OutputShape>> & phi_vals)
{
  // Start logging the dual coeff computation
  LOG_SCOPE("compute_dual_shape_coeffs()", "FE");

  const unsigned int sz=phi_vals.size();
  libmesh_error_msg_if(!sz, "ERROR: cannot compute dual shape coefficients with empty phi values");

  //compute dual basis coefficient (dual_coeff)
  dual_coeff.resize(sz, sz);
  DenseMatrix<Real> A(sz, sz), D(sz, sz);

  for (const auto i : index_range(phi_vals))
    for (const auto qp : index_range(phi_vals[i]))
    {
      D(i,i) += JxW[qp]*phi_vals[i][qp];
      for (const auto j : index_range(phi_vals))
        A(i,j) += JxW[qp]*phi_vals[i][qp]*phi_vals[j][qp];
    }

  // The integral positivity condition on the dual basis (Popp, Wohlmuth, Gee and Wall, SIAM J. Sci.
  // Comput. 34(4):B421-B446, 2012, Eq. (4.2)) fails on QUAD8 and TRI6: D(k,k) = \int N_k is -1/3 at a
  // QUAD8 corner and exactly 0 at a TRI6 vertex (their Eqs. (4.3) and (4.4)), the latter leaving those
  // dual shape functions identically zero. Biorthogonalize instead against the locally quadratic
  // transformed basis Ntilde = T N of their Sec. 4.4.1, in which each vertex absorbs a fraction alpha
  // of its adjacent mid-edge shapes. That replaces D by T^-1 diag(T d), which for this T is the sparse
  // update below. alpha = 1/5 is recommended there, and makes the weights strictly positive (QUAD8
  // 1/5 and 4/5, TRI6 1/15 and 1/10) while preserving the partition of unity. T = I elsewhere.
  // SECOND order on TRI7 is the TRI6 basis, with no degree of freedom at the centroid node, and
  // L2_LAGRANGE shares the LAGRANGE shape functions.
  if (_elem && (get_family() == LAGRANGE || get_family() == L2_LAGRANGE) &&
      (((_elem->type() == TRI6 || _elem->type() == TRI7) && sz == 6) ||
       (_elem->type() == QUAD8 && sz == 8)))
    {
      const Real alpha = Real(1)/5;
      // Mid-edge nodes are the trailing indices, and only vertex entries are written, so D(m,m) here
      // is never a value an earlier iteration modified.
      for (const auto m : make_range(_elem->n_vertices(), sz))
        for (const auto v : make_range(_elem->n_second_order_adjacent_vertices(m)))
          {
            const auto vertex = _elem->second_order_adjacent_vertex(m, v);
            D(vertex, vertex) += alpha*D(m,m);
            D(vertex, m)      -= alpha*D(m,m);
          }
    }

  // dual_coeff = A^-1*D
  for (const auto j : index_range(phi_vals))
  {
    DenseVector<Real> Dcol(sz), coeffcol(sz);
    for (const auto i : index_range(phi_vals))
      Dcol(i) = D(i, j);
    A.cholesky_solve(Dcol, coeffcol);

    for (const auto row : index_range(phi_vals))
      dual_coeff(row, j)=coeffcol(row);
  }
}

template <>
void FEGenericBase<Real>::compute_dual_shape_functions ()
{
  // Start logging the shape function computation
  LOG_SCOPE("compute_dual_shape_functions()", "FE");

  // The dual coeffs matrix should have the same size as phi
  libmesh_assert(dual_coeff.m() == phi.size());
  libmesh_assert(dual_coeff.n() == phi.size());

  // initialize dual basis
  for (const auto j : index_range(phi))
    for (const auto qp : index_range(phi[j]))
    {
      dual_phi[j][qp] = 0;
      if (calculate_dphi)
        dual_dphi[j][qp] = 0;
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2phi)
        dual_d2phi[j][qp] = 0;
#endif
    }

  // compute dual basis
  for (const auto j : index_range(phi))
    for (const auto i : index_range(phi))
      for (const auto qp : index_range(phi[j]))
      {
        dual_phi[j][qp] += dual_coeff(i, j) * phi[i][qp];
        if (calculate_dphi)
          dual_dphi[j][qp] += dual_coeff(i, j) * dphi[i][qp];
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        if (calculate_d2phi)
          dual_d2phi[j][qp] += dual_coeff(i, j) * d2phi[i][qp];
#endif
      }
}

template <typename OutputType>
void FEGenericBase<OutputType>::print_phi(std::ostream & os) const
{
  for (auto i : index_range(phi))
    for (auto j : index_range(phi[i]))
      os << " phi[" << i << "][" << j << "]=" << phi[i][j] << std::endl;
}

template <typename OutputType>
void FEGenericBase<OutputType>::print_dual_phi(std::ostream & os) const
{
  for (auto i : index_range(dual_phi))
    for (auto j : index_range(dual_phi[i]))
      os << " dual_phi[" << i << "][" << j << "]=" << dual_phi[i][j] << std::endl;
}




template <typename OutputType>
void FEGenericBase<OutputType>::print_dphi(std::ostream & os) const
{
  for (auto i : index_range(dphi))
    for (auto j : index_range(dphi[i]))
      os << " dphi[" << i << "][" << j << "]=" << dphi[i][j];
}

template <typename OutputType>
void FEGenericBase<OutputType>::print_dual_dphi(std::ostream & os) const
{
  for (auto i : index_range(dphi))
    for (auto j : index_range(dphi[i]))
      os << " dual_dphi[" << i << "][" << j << "]=" << dual_dphi[i][j];
}



template <typename OutputType>
void FEGenericBase<OutputType>::determine_calculations()
{
  this->calculations_started = true;

  // If the user did not explicitly pre-request something (or nothing)
  // to be computed, then we throw an error here.
  bool requested_ok =
    this->calculate_nothing || this->calculate_phi || this->calculate_dphi ||
    this->calculate_dphiref || this->calculate_curl_phi || this->calculate_div_phi ||
    this->calculate_map;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  requested_ok = requested_ok || this->calculate_d2phi;
#endif

  libmesh_error_msg_if(
    !requested_ok,
    "You must call one or more of the FE accessors "
    "(e.g. get_phi(), get_dphi(), get_nothing()) "
    "_before_ calling reinit()!");

  // Request whichever terms are necessary from the FEMap
  if (this->calculate_phi)
    this->_fe_trans->init_map_phi(*this);

  if (this->calculate_dphiref)
    this->_fe_trans->init_map_dphi(*this);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (this->calculate_d2phi)
    this->_fe_trans->init_map_d2phi(*this);
#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template <typename OutputType>
void FEGenericBase<OutputType>::print_d2phi(std::ostream & os) const
{
  for (auto i : index_range(dphi))
    for (auto j : index_range(dphi[i]))
      os << " d2phi[" << i << "][" << j << "]=" << d2phi[i][j];
}

template <typename OutputType>
void FEGenericBase<OutputType>::print_dual_d2phi(std::ostream & os) const
{
  for (auto i : index_range(dual_d2phi))
    for (auto j : index_range(dual_d2phi[i]))
      os << " dual_d2phi[" << i << "][" << j << "]=" << dual_d2phi[i][j];
}

#endif



#ifdef LIBMESH_ENABLE_AMR

template <typename OutputType>
void
FEGenericBase<OutputType>::coarsened_dof_values(const NumericVector<Number> & old_vector,
                                                const DofMap & dof_map,
                                                const Elem * elem,
                                                DenseVector<Number> & Ue,
                                                const unsigned int var,
                                                const bool use_old_dof_indices)
{
  // Side/edge local DOF indices
  std::vector<unsigned int> new_side_dofs, old_side_dofs;

  // FIXME: what about 2D shells in 3D space?
  unsigned int dim = elem->dim();

  // Cache n_children(); it's a virtual call but it's const.
  const unsigned int n_children = elem->n_children();

  // We use local FE objects for now
  // FIXME: we should use more, external objects instead for efficiency
  const FEType & base_fe_type = dof_map.variable_type(var);
  std::unique_ptr<FEGenericBase<OutputShape>> fe
    (FEGenericBase<OutputShape>::build(dim, base_fe_type));
  std::unique_ptr<FEGenericBase<OutputShape>> fe_coarse
    (FEGenericBase<OutputShape>::build(dim, base_fe_type));

  std::unique_ptr<QBase> qrule     (base_fe_type.default_quadrature_rule(dim));
  std::unique_ptr<QBase> qedgerule (base_fe_type.default_quadrature_rule(1));
  std::unique_ptr<QBase> qsiderule (base_fe_type.default_quadrature_rule(dim-1));
  std::vector<Point> coarse_qpoints;

  // The values of the shape functions at the quadrature
  // points
  const std::vector<std::vector<OutputShape>> & phi_values =
    fe->get_phi();
  const std::vector<std::vector<OutputShape>> & phi_coarse =
    fe_coarse->get_phi();

  // The gradients of the shape functions at the quadrature
  // points on the child element.
  const std::vector<std::vector<OutputGradient>> * dphi_values =
    nullptr;
  const std::vector<std::vector<OutputGradient>> * dphi_coarse =
    nullptr;

  const FEContinuity cont = fe->get_continuity();

  if (cont == C_ONE)
    {
      const std::vector<std::vector<OutputGradient>> &
        ref_dphi_values = fe->get_dphi();
      dphi_values = &ref_dphi_values;
      const std::vector<std::vector<OutputGradient>> &
        ref_dphi_coarse = fe_coarse->get_dphi();
      dphi_coarse = &ref_dphi_coarse;
    }

  // The Jacobian * quadrature weight at the quadrature points
  const std::vector<Real> & JxW =
    fe->get_JxW();

  // The XYZ locations of the quadrature points on the
  // child element
  const std::vector<Point> & xyz_values =
    fe->get_xyz();

  // Number of nodes on parent element
  const unsigned int n_nodes = elem->n_nodes();

  // Number of dofs on parent element
  const unsigned int new_n_dofs =
    FEInterface::n_dofs(base_fe_type, elem->max_descendant_p_level(), elem);

  // Fixed vs. free DoFs on edge/face projections
  std::vector<char> dof_is_fixed(new_n_dofs, false); // bools
  std::vector<int> free_dof(new_n_dofs, 0);

  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  Ue.resize(new_n_dofs); Ue.zero();


  // When coarsening, in general, we need a series of
  // projections to ensure a unique and continuous
  // solution.  We start by interpolating nodes, then
  // hold those fixed and project edges, then
  // hold those fixed and project faces, then
  // hold those fixed and project interiors

  // Copy node values first
  {
    std::vector<dof_id_type> node_dof_indices;
    if (use_old_dof_indices)
      dof_map.old_dof_indices (elem, node_dof_indices, var);
    else
      dof_map.dof_indices (elem, node_dof_indices, var);

    unsigned int current_dof = 0;
    for (unsigned int n=0; n!= n_nodes; ++n)
      {
        // FIXME: this should go through the DofMap,
        // not duplicate dof_indices code badly!
        const unsigned int my_nc =
          FEInterface::n_dofs_at_node (base_fe_type, elem->max_descendant_p_level(), elem, n);
        if (!elem->is_vertex(n))
          {
            current_dof += my_nc;
            continue;
          }

        // We're assuming here that child n shares vertex n,
        // which is wrong on non-simplices right now
        // ... but this code isn't necessary except on elements
        // where p refinement creates more vertex dofs; we have
        // no such elements yet.
        int extra_order = 0;
        // if (elem->child_ptr(n)->p_level() < elem->p_level())
        //   extra_order = elem->child_ptr(n)->p_level();
        const unsigned int nc =
          FEInterface::n_dofs_at_node (base_fe_type, extra_order, elem, n);
        for (unsigned int i=0; i!= nc; ++i)
          {
            Ue(current_dof) =
              old_vector(node_dof_indices[current_dof]);
            dof_is_fixed[current_dof] = true;
            current_dof++;
          }
      }
  }

  FEType fe_type = base_fe_type, temp_fe_type;
  fe_type.order = fe_type.order + elem->max_descendant_p_level();

  // In 3D, project any edge values next
  if (dim > 2 && cont != DISCONTINUOUS)
    for (auto e : elem->edge_index_range())
      {
        FEInterface::dofs_on_edge(elem, dim, fe_type,
                                  e, new_side_dofs);

        const unsigned int n_new_side_dofs =
          cast_int<unsigned int>(new_side_dofs.size());

        // Some edge dofs are on nodes and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (unsigned int i=0; i != n_new_side_dofs; ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new edge coefficients
        DenseVector<Number> Uedge(free_dofs);

        // Add projection terms from each child sharing
        // this edge
        for (unsigned int c=0; c != n_children; ++c)
          {
            if (!elem->is_child_on_edge(c,e))
              continue;
            const Elem * child = elem->child_ptr(c);

            std::vector<dof_id_type> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                                       child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                                   child_dof_indices, var);
            const unsigned int child_n_dofs =
              cast_int<unsigned int>
              (child_dof_indices.size());

            temp_fe_type = base_fe_type;
            temp_fe_type.order = temp_fe_type.order + child->p_level();

            FEInterface::dofs_on_edge(child, dim,
                                      temp_fe_type, e, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's edge
            fe->attach_quadrature_rule (qedgerule.get());
            fe->edge_reinit (child, e);
            const unsigned int n_qp = qedgerule->n_points();

            FEMap::inverse_map (dim, elem, xyz_values,
                                coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                OutputNumber fineval = libMesh::zero;
                // solution grad at the quadrature point
                OutputNumberGradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      (old_vector(child_dof_indices[i])*
                       phi_values[i][qp]);
                    if (cont == C_ONE)
                      finegrad += (*dphi_values)[i][qp] *
                        old_vector(child_dof_indices[i]);
                  }

                // Form edge projection matrix
                for (unsigned int sidei=0, freei=0; sidei != n_new_side_dofs; ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0; sidej != n_new_side_dofs; ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) +=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) +=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += TensorTools::inner_product(phi_coarse[i][qp],
                                                            fineval) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uedge);

        // Transfer new edge solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number & ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uedge(i)) < TOLERANCE);
            ui = Uedge(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] = true;
          }
      }

  // Project any side values (edges in 2D, faces in 3D)
  if (dim > 1 && cont != DISCONTINUOUS)
    for (auto s : elem->side_index_range())
      {
        FEInterface::dofs_on_side(elem, dim, fe_type,
                                  s, new_side_dofs);

        const unsigned int n_new_side_dofs =
          cast_int<unsigned int>(new_side_dofs.size());

        // Some side dofs are on nodes/edges and already
        // fixed, others are free to calculate
        unsigned int free_dofs = 0;
        for (unsigned int i=0; i != n_new_side_dofs; ++i)
          if (!dof_is_fixed[new_side_dofs[i]])
            free_dof[free_dofs++] = i;
        Ke.resize (free_dofs, free_dofs); Ke.zero();
        Fe.resize (free_dofs); Fe.zero();
        // The new side coefficients
        DenseVector<Number> Uside(free_dofs);

        // Add projection terms from each child sharing
        // this side
        for (unsigned int c=0; c != n_children; ++c)
          {
            if (!elem->is_child_on_side(c,s))
              continue;
            const Elem * child = elem->child_ptr(c);

            std::vector<dof_id_type> child_dof_indices;
            if (use_old_dof_indices)
              dof_map.old_dof_indices (child,
                                       child_dof_indices, var);
            else
              dof_map.dof_indices (child,
                                   child_dof_indices, var);
            const unsigned int child_n_dofs =
              cast_int<unsigned int>
              (child_dof_indices.size());

            temp_fe_type = base_fe_type;
            temp_fe_type.order = temp_fe_type.order + child->p_level();

            FEInterface::dofs_on_side(child, dim,
                                      temp_fe_type, s, old_side_dofs);

            // Initialize both child and parent FE data
            // on the child's side
            fe->attach_quadrature_rule (qsiderule.get());
            fe->reinit (child, s);
            const unsigned int n_qp = qsiderule->n_points();

            FEMap::inverse_map (dim, elem, xyz_values,
                                coarse_qpoints);

            fe_coarse->reinit(elem, &coarse_qpoints);

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution value at the quadrature point
                OutputNumber fineval = libMesh::zero;
                // solution grad at the quadrature point
                OutputNumberGradient finegrad;

                // Sum the solution values * the DOF
                // values at the quadrature point to
                // get the solution value and gradient.
                for (unsigned int i=0; i<child_n_dofs;
                     i++)
                  {
                    fineval +=
                      old_vector(child_dof_indices[i]) *
                      phi_values[i][qp];
                    if (cont == C_ONE)
                      finegrad += (*dphi_values)[i][qp] *
                        old_vector(child_dof_indices[i]);
                  }

                // Form side projection matrix
                for (unsigned int sidei=0, freei=0; sidei != n_new_side_dofs; ++sidei)
                  {
                    unsigned int i = new_side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0; sidej != n_new_side_dofs; ++sidej)
                      {
                        unsigned int j =
                          new_side_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) +=
                            TensorTools::inner_product(phi_coarse[i][qp],
                                                       phi_coarse[j][qp]) *
                            JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) +=
                                TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                           (*dphi_coarse)[j][qp]) *
                                JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += TensorTools::inner_product(fineval, phi_coarse[i][qp]) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) +=
                        TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
                    freei++;
                  }
              }
          }
        Ke.cholesky_solve(Fe, Uside);

        // Transfer new side solutions to element
        for (unsigned int i=0; i != free_dofs; ++i)
          {
            Number & ui = Ue(new_side_dofs[free_dof[i]]);
            libmesh_assert(std::abs(ui) < TOLERANCE ||
                           std::abs(ui - Uside(i)) < TOLERANCE);
            ui = Uside(i);
            dof_is_fixed[new_side_dofs[free_dof[i]]] = true;
          }
      }

  // Project the interior values, finally

  // Some interior dofs are on nodes/edges/sides and
  // already fixed, others are free to calculate
  unsigned int free_dofs = 0;
  for (unsigned int i=0; i != new_n_dofs; ++i)
    if (!dof_is_fixed[i])
      free_dof[free_dofs++] = i;
  Ke.resize (free_dofs, free_dofs); Ke.zero();
  Fe.resize (free_dofs); Fe.zero();
  // The new interior coefficients
  DenseVector<Number> Uint(free_dofs);

  // Add projection terms from each child
  for (auto & child : elem->child_ref_range())
    {
      std::vector<dof_id_type> child_dof_indices;
      if (use_old_dof_indices)
        dof_map.old_dof_indices (&child,
                                 child_dof_indices, var);
      else
        dof_map.dof_indices (&child,
                             child_dof_indices, var);
      const unsigned int child_n_dofs =
        cast_int<unsigned int>
        (child_dof_indices.size());

      // Initialize both child and parent FE data
      // on the child's quadrature points
      fe->attach_quadrature_rule (qrule.get());
      fe->reinit (&child);
      const unsigned int n_qp = qrule->n_points();

      FEMap::inverse_map (dim, elem, xyz_values, coarse_qpoints);

      fe_coarse->reinit(elem, &coarse_qpoints);

      // Loop over the quadrature points
      for (unsigned int qp=0; qp<n_qp; qp++)
        {
          // solution value at the quadrature point
          OutputNumber fineval = libMesh::zero;
          // solution grad at the quadrature point
          OutputNumberGradient finegrad;

          // Sum the solution values * the DOF
          // values at the quadrature point to
          // get the solution value and gradient.
          for (unsigned int i=0; i<child_n_dofs; i++)
            {
              fineval +=
                (old_vector(child_dof_indices[i]) *
                 phi_values[i][qp]);
              if (cont == C_ONE)
                finegrad += (*dphi_values)[i][qp] *
                  old_vector(child_dof_indices[i]);
            }

          // Form interior projection matrix
          for (unsigned int i=0, freei=0;
               i != new_n_dofs; ++i)
            {
              // fixed DoFs aren't test functions
              if (dof_is_fixed[i])
                continue;
              for (unsigned int j=0, freej=0; j !=
                     new_n_dofs; ++j)
                {
                  if (dof_is_fixed[j])
                    Fe(freei) -=
                      TensorTools::inner_product(phi_coarse[i][qp],
                                                 phi_coarse[j][qp]) *
                      JxW[qp] * Ue(j);
                  else
                    Ke(freei,freej) +=
                      TensorTools::inner_product(phi_coarse[i][qp],
                                                 phi_coarse[j][qp]) *
                      JxW[qp];
                  if (cont == C_ONE)
                    {
                      if (dof_is_fixed[j])
                        Fe(freei) -=
                          TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                     (*dphi_coarse)[j][qp]) *
                          JxW[qp] * Ue(j);
                      else
                        Ke(freei,freej) +=
                          TensorTools::inner_product((*dphi_coarse)[i][qp],
                                                     (*dphi_coarse)[j][qp]) *
                          JxW[qp];
                    }
                  if (!dof_is_fixed[j])
                    freej++;
                }
              Fe(freei) += TensorTools::inner_product(phi_coarse[i][qp], fineval) *
                JxW[qp];
              if (cont == C_ONE)
                Fe(freei) += TensorTools::inner_product(finegrad, (*dphi_coarse)[i][qp]) * JxW[qp];
              freei++;
            }
        }
    }
  Ke.cholesky_solve(Fe, Uint);

  // Transfer new interior solutions to element
  for (unsigned int i=0; i != free_dofs; ++i)
    {
      Number & ui = Ue(free_dof[i]);
      libmesh_assert(std::abs(ui) < TOLERANCE ||
                     std::abs(ui - Uint(i)) < TOLERANCE);
      ui = Uint(i);
      // We should be fixing all dofs by now; no need to keep track of
      // that unless we're debugging
#ifndef NDEBUG
      dof_is_fixed[free_dof[i]] = true;
#endif
    }

#ifndef NDEBUG
  // Make sure every DoF got reached!
  for (unsigned int i=0; i != new_n_dofs; ++i)
    libmesh_assert(dof_is_fixed[i]);
#endif
}



template <typename OutputType>
void
FEGenericBase<OutputType>::coarsened_dof_values(const NumericVector<Number> & old_vector,
                                                const DofMap & dof_map,
                                                const Elem * elem,
                                                DenseVector<Number> & Ue,
                                                const bool use_old_dof_indices)
{
  Ue.resize(0);

  for (auto v : make_range(dof_map.n_variables()))
    {
      DenseVector<Number> Usub;

      coarsened_dof_values(old_vector, dof_map, elem, Usub,
                           v, use_old_dof_indices);

      Ue.append (Usub);
    }
}



template <typename OutputType>
void
FEGenericBase<OutputType>::compute_proj_constraints (DofConstraints & constraints,
                                                     DofMap & dof_map,
                                                     const unsigned int variable_number,
                                                     const Elem * elem)
{
  libmesh_assert(elem);

  const unsigned int Dim = elem->dim();

  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  const Variable & var = dof_map.variable(variable_number);
  const FEType & base_fe_type = var.type();
  const bool add_p_level = base_fe_type.p_refinement;

  // Construct FE objects for this element and its neighbors.
  std::unique_ptr<FEGenericBase<OutputShape>> my_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));
  my_fe->add_p_level_in_reinit(add_p_level);
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE ||
                  cont == SIDE_DISCONTINUOUS);

  // this would require some generalisation:
  //  - e.g. the 'my_fe'-object needs generalisation
  //  - due to lack of one-to-one correspondence of DOFs and nodes,
  //    this doesn't work easily.
  if (elem->infinite())
    libmesh_not_implemented();

  std::unique_ptr<FEGenericBase<OutputShape>> neigh_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));
  neigh_fe->add_p_level_in_reinit(add_p_level);

  QGauss my_qface(Dim-1, base_fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> neigh_qface;

  const std::vector<Real> & JxW = my_fe->get_JxW();
  const std::vector<Point> & q_point = my_fe->get_xyz();
  const std::vector<std::vector<OutputShape>> & phi = my_fe->get_phi();
  const std::vector<std::vector<OutputShape>> & neigh_phi =
    neigh_fe->get_phi();
  const std::vector<Point> * face_normals = nullptr;
  const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
  const std::vector<std::vector<OutputGradient>> * neigh_dphi = nullptr;

  std::vector<dof_id_type> my_dof_indices, neigh_dof_indices;
  std::vector<unsigned int> my_side_dofs, neigh_side_dofs;

  if (cont == C_ONE)
    {
      const std::vector<Point> & ref_face_normals =
        my_fe->get_normals();
      face_normals = &ref_face_normals;
      const std::vector<std::vector<OutputGradient>> & ref_dphi =
        my_fe->get_dphi();
      dphi = &ref_dphi;
      const std::vector<std::vector<OutputGradient>> & ref_neigh_dphi =
        neigh_fe->get_dphi();
      neigh_dphi = &ref_neigh_dphi;
    }

  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;
  std::vector<DenseVector<Real>> Ue;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (auto s : elem->side_index_range())
    {
      // Get pointers to the element's neighbor.
      const Elem * neigh = elem->neighbor_ptr(s);

      if (!neigh)
        continue;

      if (!var.active_on_subdomain(neigh->subdomain_id()))
        continue;

      // h refinement constraints:
      // constrain dofs shared between
      // this element and ones coarser
      // than this element.
      if (neigh->level() < elem->level())
        {
          unsigned int s_neigh = neigh->which_neighbor_am_i(elem);
          libmesh_assert_less (s_neigh, neigh->n_neighbors());

          // Find the minimum p level; we build the h constraint
          // matrix with this and then constrain away all higher p
          // DoFs.
          libmesh_assert(neigh->active());
          const unsigned int min_p_level = add_p_level *
            std::min(elem->p_level(), neigh->p_level());
          // we may need to make the FE objects reinit with the
          // minimum shared p_level
          const unsigned int old_elem_level = add_p_level * elem->p_level();
          if (old_elem_level != min_p_level)
            my_fe->set_fe_order(my_fe->get_fe_type().order.get_order() + min_p_level - old_elem_level);
          const unsigned int old_neigh_level = add_p_level * neigh->p_level();
          if (old_neigh_level != min_p_level)
            neigh_fe->set_fe_order(neigh_fe->get_fe_type().order.get_order() + min_p_level - old_neigh_level);

          my_fe->reinit(elem, s);

          // This function gets called element-by-element, so there
          // will be a lot of memory allocation going on.  We can
          // at least minimize this for the case of the dof indices
          // by efficiently preallocating the requisite storage.
          // n_nodes is not necessarily n_dofs, but it is better
          // than nothing!
          my_dof_indices.reserve    (elem->n_nodes());
          neigh_dof_indices.reserve (neigh->n_nodes());

          dof_map.dof_indices (elem, my_dof_indices,
                               variable_number,
                               min_p_level);
          dof_map.dof_indices (neigh, neigh_dof_indices,
                               variable_number,
                               min_p_level);

          const unsigned int n_qp = my_qface.n_points();

          FEMap::inverse_map (Dim, neigh, q_point, neigh_qface);

          neigh_fe->reinit(neigh, &neigh_qface);

          // We're only concerned with DOFs whose values (and/or first
          // derivatives for C1 elements) are supported on side nodes
          FEType elem_fe_type = base_fe_type;
          if (old_elem_level != min_p_level)
            elem_fe_type.order = base_fe_type.order.get_order() + min_p_level - old_elem_level;
          FEType neigh_fe_type = base_fe_type;
          if (old_neigh_level != min_p_level)
            neigh_fe_type.order = base_fe_type.order.get_order() + min_p_level - old_neigh_level;
          FEInterface::dofs_on_side(elem,  Dim, elem_fe_type,  s,       my_side_dofs);
          FEInterface::dofs_on_side(neigh, Dim, neigh_fe_type, s_neigh, neigh_side_dofs);

          const unsigned int n_side_dofs =
            cast_int<unsigned int>(my_side_dofs.size());
          libmesh_assert_equal_to (n_side_dofs, neigh_side_dofs.size());

#ifndef NDEBUG
          for (auto i : my_side_dofs)
            libmesh_assert_less(i, my_dof_indices.size());
          for (auto i : neigh_side_dofs)
            libmesh_assert_less(i, neigh_dof_indices.size());
#endif

          Ke.resize (n_side_dofs, n_side_dofs);
          Ue.resize(n_side_dofs);

          // Form the projection matrix, (inner product of fine basis
          // functions against fine test functions)
          for (unsigned int is = 0; is != n_side_dofs; ++is)
            {
              const unsigned int i = my_side_dofs[is];
              for (unsigned int js = 0; js != n_side_dofs; ++js)
                {
                  const unsigned int j = my_side_dofs[js];
                  for (unsigned int qp = 0; qp != n_qp; ++qp)
                    {
                      Ke(is,js) += JxW[qp] * TensorTools::inner_product(phi[i][qp], phi[j][qp]);
                      if (cont == C_ONE)
                        Ke(is,js) += JxW[qp] *
                          TensorTools::inner_product((*dphi)[i][qp] *
                                                     (*face_normals)[qp],
                                                     (*dphi)[j][qp] *
                                                     (*face_normals)[qp]);
                    }
                }
            }

          // Form the right hand sides, (inner product of coarse basis
          // functions against fine test functions)
          for (unsigned int is = 0; is != n_side_dofs; ++is)
            {
              const unsigned int i = neigh_side_dofs[is];
              Fe.resize (n_side_dofs);
              for (unsigned int js = 0; js != n_side_dofs; ++js)
                {
                  const unsigned int j = my_side_dofs[js];
                  for (unsigned int qp = 0; qp != n_qp; ++qp)
                    {
                      Fe(js) += JxW[qp] *
                        TensorTools::inner_product(neigh_phi[i][qp],
                                                   phi[j][qp]);
                      if (cont == C_ONE)
                        Fe(js) += JxW[qp] *
                          TensorTools::inner_product((*neigh_dphi)[i][qp] *
                                                     (*face_normals)[qp],
                                                     (*dphi)[j][qp] *
                                                     (*face_normals)[qp]);
                    }
                }
              Ke.cholesky_solve(Fe, Ue[is]);
            }

          for (unsigned int js = 0; js != n_side_dofs; ++js)
            {
              const unsigned int j = my_side_dofs[js];
              const dof_id_type my_dof_g = my_dof_indices[j];
              libmesh_assert_not_equal_to (my_dof_g, DofObject::invalid_id);

              // Hunt for "constraining against myself" cases before
              // we bother creating a constraint row
              bool self_constraint = false;
              for (unsigned int is = 0; is != n_side_dofs; ++is)
                {
                  const unsigned int i = neigh_side_dofs[is];
                  const dof_id_type their_dof_g = neigh_dof_indices[i];
                  libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);

                  if (their_dof_g == my_dof_g)
                    {
#ifndef NDEBUG
                      const Real their_dof_value = Ue[is](js);
                      libmesh_assert_less (std::abs(their_dof_value-1.),
                                           10*TOLERANCE);

                      for (unsigned int k = 0; k != n_side_dofs; ++k)
                        libmesh_assert(k == is ||
                                       std::abs(Ue[k](js)) <
                                       10*TOLERANCE);
#endif

                      self_constraint = true;
                      break;
                    }
                }

              if (self_constraint)
                continue;

              DofConstraintRow * constraint_row;

              // we may be running constraint methods concurrently
              // on multiple threads, so we need a lock to
              // ensure that this constraint is "ours"
              {
                Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                if (dof_map.is_constrained_dof(my_dof_g))
                  continue;

                constraint_row = &(constraints[my_dof_g]);
                libmesh_assert(constraint_row->empty());
              }

              for (unsigned int is = 0; is != n_side_dofs; ++is)
                {
                  const unsigned int i = neigh_side_dofs[is];
                  const dof_id_type their_dof_g = neigh_dof_indices[i];
                  libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);
                  libmesh_assert_not_equal_to (their_dof_g, my_dof_g);

                  const Real their_dof_value = Ue[is](js);

                  if (std::abs(their_dof_value) < 10*TOLERANCE)
                    continue;

                  constraint_row->emplace(their_dof_g, their_dof_value);
                }
            }

          my_fe->set_fe_order(my_fe->get_fe_type().order.get_order() + old_elem_level - min_p_level);
          neigh_fe->set_fe_order(neigh_fe->get_fe_type().order.get_order() + old_neigh_level - min_p_level);
        }

      if (add_p_level)
      {
        // p refinement constraints:
        // constrain dofs shared between
        // active elements and neighbors with
        // lower polynomial degrees
        const unsigned int min_p_level =
          neigh->min_p_level_by_neighbor(elem, elem->p_level());
        if (min_p_level < elem->p_level())
        {
          // Adaptive p refinement of non-hierarchic bases will
          // require more coding
          libmesh_assert(my_fe->is_hierarchic());
          dof_map.constrain_p_dofs(variable_number, elem,
                                   s, min_p_level);
        }
      }
    }
}

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC
template <typename OutputType>
void
FEGenericBase<OutputType>::
compute_periodic_constraints (DofConstraints & constraints,
                              DofMap & dof_map,
                              const PeriodicBoundaries & boundaries,
                              const MeshBase & mesh,
                              const PointLocatorBase * point_locator,
                              const unsigned int variable_number,
                              const Elem * elem)
{
  // Only bother if we truly have periodic boundaries
  if (boundaries.empty())
    return;

  libmesh_assert(elem);

  // Only constrain active elements with this method
  if (!elem->active())
    return;

  if (elem->infinite())
    libmesh_not_implemented();

  const unsigned int Dim = elem->dim();

  // We need sys_number and variable_number for DofObject methods
  // later
  const unsigned int sys_number = dof_map.sys_number();

  const FEType & base_fe_type = dof_map.variable_type(variable_number);

  // Construct FE objects for this element and its pseudo-neighbors.
  std::unique_ptr<FEGenericBase<OutputShape>> my_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));
  const FEContinuity cont = my_fe->get_continuity();

  // We don't need to constrain discontinuous elements
  if (cont == DISCONTINUOUS)
    return;
  libmesh_assert (cont == C_ZERO || cont == C_ONE);

  // We'll use element size to generate relative tolerances later
  const Real elem_hmin = elem->hmin();

  std::unique_ptr<FEGenericBase<OutputShape>> neigh_fe
    (FEGenericBase<OutputShape>::build(Dim, base_fe_type));

  QGauss my_qface(Dim-1, base_fe_type.default_quadrature_order());
  my_fe->attach_quadrature_rule (&my_qface);
  std::vector<Point> neigh_qface;

  const std::vector<Real> & JxW = my_fe->get_JxW();
  const std::vector<Point> & q_point = my_fe->get_xyz();
  const std::vector<std::vector<OutputShape>> & phi = my_fe->get_phi();
  const std::vector<std::vector<OutputShape>> & neigh_phi =
    neigh_fe->get_phi();
  const std::vector<Point> * face_normals = nullptr;
  const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
  const std::vector<std::vector<OutputGradient>> * neigh_dphi = nullptr;
  std::vector<dof_id_type> my_dof_indices, neigh_dof_indices;
  std::vector<unsigned int> my_side_dofs, neigh_side_dofs;

  if (cont != C_ZERO)
    {
      const std::vector<Point> & ref_face_normals =
        my_fe->get_normals();
      face_normals = &ref_face_normals;
      const std::vector<std::vector<OutputGradient>> & ref_dphi =
        my_fe->get_dphi();
      dphi = &ref_dphi;
      const std::vector<std::vector<OutputGradient>> & ref_neigh_dphi =
        neigh_fe->get_dphi();
      neigh_dphi = &ref_neigh_dphi;
    }

  DenseMatrix<Real> Ke;
  DenseVector<Real> Fe;
  std::vector<DenseVector<Real>> Ue;

  // Container to catch the boundary ids that BoundaryInfo hands us.
  std::vector<boundary_id_type> bc_ids;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  const unsigned short int max_ns = elem->n_sides();
  for (unsigned short int s = 0; s != max_ns; ++s)
    {
      if (elem->neighbor_ptr(s))
        continue;

      mesh.get_boundary_info().boundary_ids (elem, s, bc_ids);

      for (const auto & boundary_id : bc_ids)
        {
          const PeriodicBoundaryBase * periodic = boundaries.boundary(boundary_id);
          if (!periodic || !periodic->is_my_variable(variable_number))
            continue;

          libmesh_assert(point_locator);

          // Get pointers to the element's neighbor.
          unsigned int s_neigh;
          const Elem * neigh = boundaries.neighbor(boundary_id, *point_locator, elem, s, &s_neigh);

          libmesh_error_msg_if(neigh == nullptr,
                               "PeriodicBoundaries point locator object returned nullptr!");

          // periodic (and possibly h refinement) constraints:
          // constrain dofs shared between
          // this element and ones as coarse
          // as or coarser than this element.
          if (neigh->level() <= elem->level())
            {
#ifdef LIBMESH_ENABLE_AMR
              // Find the minimum p level; we build the h constraint
              // matrix with this and then constrain away all higher p
              // DoFs.
              libmesh_assert(neigh->active());
              const unsigned int min_p_level =
                std::min(elem->p_level(), neigh->p_level());

              // we may need to make the FE objects reinit with the
              // minimum shared p_level
              // FIXME - I hate using const_cast<> and avoiding
              // accessor functions; there's got to be a
              // better way to do this!
              const unsigned int old_elem_level = elem->p_level();
              if (old_elem_level != min_p_level)
                (const_cast<Elem *>(elem))->hack_p_level(min_p_level);
              const unsigned int old_neigh_level = neigh->p_level();
              if (old_neigh_level != min_p_level)
                (const_cast<Elem *>(neigh))->hack_p_level(min_p_level);
#endif // #ifdef LIBMESH_ENABLE_AMR

              // We can do a projection with a single integration,
              // due to the assumption of nested finite element
              // subspaces.
              // FIXME: it might be more efficient to do nodes,
              // then edges, then side, to reduce the size of the
              // Cholesky factorization(s)
              my_fe->reinit(elem, s);

              dof_map.dof_indices (elem, my_dof_indices,
                                   variable_number);
              dof_map.dof_indices (neigh, neigh_dof_indices,
                                   variable_number);

              // We use neigh_dof_indices_all_variables in the case that the
              // periodic boundary condition involves mappings between multiple
              // variables.
              std::vector<std::vector<dof_id_type>> neigh_dof_indices_all_variables;
              if(periodic->has_transformation_matrix())
                {
                  const std::set<unsigned int> & variables = periodic->get_variables();
                  neigh_dof_indices_all_variables.resize(variables.size());
                  unsigned int index = 0;
                  for(unsigned int var : variables)
                    {
                      dof_map.dof_indices (neigh, neigh_dof_indices_all_variables[index],
                                           var);
                      index++;
                    }
                }

              const unsigned int n_qp = my_qface.n_points();

              // Translate the quadrature points over to the
              // neighbor's boundary
              std::vector<Point> neigh_point(q_point.size());
              for (auto i : index_range(neigh_point))
                neigh_point[i] = periodic->get_corresponding_pos(q_point[i]);

              FEMap::inverse_map (Dim, neigh, neigh_point,
                                  neigh_qface);

              neigh_fe->reinit(neigh, &neigh_qface);

              // We're only concerned with DOFs whose values (and/or first
              // derivatives for C1 elements) are supported on side nodes
              FEInterface::dofs_on_side(elem, Dim, base_fe_type, s, my_side_dofs);
              FEInterface::dofs_on_side(neigh, Dim, base_fe_type, s_neigh, neigh_side_dofs);

              // We're done with functions that examine Elem::p_level(),
              // so let's unhack those levels
#ifdef LIBMESH_ENABLE_AMR
              if (elem->p_level() != old_elem_level)
                (const_cast<Elem *>(elem))->hack_p_level(old_elem_level);
              if (neigh->p_level() != old_neigh_level)
                (const_cast<Elem *>(neigh))->hack_p_level(old_neigh_level);
#endif // #ifdef LIBMESH_ENABLE_AMR

              const unsigned int n_side_dofs =
                cast_int<unsigned int>
                (my_side_dofs.size());
              libmesh_assert_equal_to (n_side_dofs, neigh_side_dofs.size());

              Ke.resize (n_side_dofs, n_side_dofs);
              Ue.resize(n_side_dofs);

              // Form the projection matrix, (inner product of fine basis
              // functions against fine test functions)
              for (unsigned int is = 0; is != n_side_dofs; ++is)
                {
                  const unsigned int i = my_side_dofs[is];
                  for (unsigned int js = 0; js != n_side_dofs; ++js)
                    {
                      const unsigned int j = my_side_dofs[js];
                      for (unsigned int qp = 0; qp != n_qp; ++qp)
                        {
                          Ke(is,js) += JxW[qp] *
                            TensorTools::inner_product(phi[i][qp],
                                                       phi[j][qp]);
                          if (cont != C_ZERO)
                            Ke(is,js) += JxW[qp] *
                              TensorTools::inner_product((*dphi)[i][qp] *
                                                         (*face_normals)[qp],
                                                         (*dphi)[j][qp] *
                                                         (*face_normals)[qp]);
                        }
                    }
                }

              // Form the right hand sides, (inner product of coarse basis
              // functions against fine test functions)
              for (unsigned int is = 0; is != n_side_dofs; ++is)
                {
                  const unsigned int i = neigh_side_dofs[is];
                  Fe.resize (n_side_dofs);
                  for (unsigned int js = 0; js != n_side_dofs; ++js)
                    {
                      const unsigned int j = my_side_dofs[js];
                      for (unsigned int qp = 0; qp != n_qp; ++qp)
                        {
                          Fe(js) += JxW[qp] *
                            TensorTools::inner_product(neigh_phi[i][qp],
                                                       phi[j][qp]);
                          if (cont != C_ZERO)
                            Fe(js) += JxW[qp] *
                              TensorTools::inner_product((*neigh_dphi)[i][qp] *
                                                         (*face_normals)[qp],
                                                         (*dphi)[j][qp] *
                                                         (*face_normals)[qp]);
                        }
                    }
                  Ke.cholesky_solve(Fe, Ue[is]);
                }

              // Make sure we're not adding cyclic constraints due to
              // the redundancy in the way we add periodic boundary
              // constraints
              //
              // In order for this to work while threaded or on
              // distributed meshes, we need a rigorous way to
              // avoid cycles of constraints.  On distributed
              // meshes with remote point-neighbors in particular, not
              // every rank may be able to see the "root"
              // unconstrained DoF(s) that should constrain its own
              // periodic constrained DoFs.
              //
              // We handle this via a total ordering on nodes (which
              // is where all boundary DoFs that might have a periodic
              // boundary live).
              //
              // If the paired element we're looking at is coarser
              // than our element, we can constrain to it.
              //
              // If the paired active neighbor we're looking at is at
              // the same level as our active element, then it is
              // "conforming" across the periodic boundary, and our
              // boundary nodes should match theirs one-to-one; we
              // will just use nodes' lexicographic ordering to
              // determine which constrains which.
              std::set<dof_id_type> dofs_to_not_constrain;

              if (neigh->level() == elem->level())
                {
                  auto elem_nodes = elem->nodes_on_side(s);
                  auto neigh_nodes = neigh->nodes_on_side(s_neigh);

                  const Real tol = elem_hmin * TOLERANCE;

                  auto do_not_constrain =
                    [&dofs_to_not_constrain, sys_number, variable_number]
                    (const DofObject & d)
                    {
                      for (auto i : make_range(d.n_comp(sys_number,
                                                        variable_number)))
                        dofs_to_not_constrain.insert
                          (d.dof_number (sys_number, variable_number, i));
                    };

                  for (auto n : elem_nodes)
                    {
                      const Node & my_node = elem->node_ref(n);

                      const Point neigh_pt =
                        periodic->get_corresponding_pos(my_node);

                      if (my_node.absolute_fuzzy_equals(neigh_pt, tol))
                        {
                          do_not_constrain(my_node);
                          continue;
                        }

                      for (auto nn : neigh_nodes)
                        {
                          const Point & neigh_node = neigh->point(nn);

                          if (neigh_node.absolute_fuzzy_equals(neigh_pt, tol))
                            {
                              if (neigh_node > my_node)
                                do_not_constrain(my_node);
                              break;
                            }
                        }
                    }
                }

              for (unsigned int js = 0; js != n_side_dofs; ++js)
                {
                  const unsigned int j = my_side_dofs[js];
                  const dof_id_type my_dof_g = my_dof_indices[j];
                  libmesh_assert_not_equal_to (my_dof_g, DofObject::invalid_id);

                  if (dofs_to_not_constrain.count(my_dof_g))
                    continue;

                  DofConstraintRow * constraint_row;

                  // we may be running constraint methods concurrently
                  // on multiple threads, so we need a lock to
                  // ensure that this constraint is "ours"
                  {
                    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                    if (dof_map.is_constrained_dof(my_dof_g))
                      continue;

                    constraint_row = &(constraints[my_dof_g]);
                    libmesh_assert(constraint_row->empty());
                  }

                  for (unsigned int is = 0; is != n_side_dofs; ++is)
                    {
                      const unsigned int i = neigh_side_dofs[is];
                      const dof_id_type their_dof_g = neigh_dof_indices[i];
                      libmesh_assert_not_equal_to (their_dof_g, DofObject::invalid_id);

                      // Periodic constraints should never be
                      // self-constraints
                      // libmesh_assert_not_equal_to (their_dof_g, my_dof_g);

                      const Real their_dof_value = Ue[is](js);

                      if (their_dof_g == my_dof_g)
                        {
                          libmesh_assert_less (std::abs(their_dof_value-1.), 1.e-5);
                          for (unsigned int k = 0; k != n_side_dofs; ++k)
                            libmesh_assert(k == is || std::abs(Ue[k](js)) < 1.e-5);
                          continue;
                        }

                      if (std::abs(their_dof_value) < 10*TOLERANCE)
                        continue;

                      if(!periodic->has_transformation_matrix())
                        {
                          constraint_row->emplace(their_dof_g, their_dof_value);
                        }
                      else
                        {
                          // In this case the current variable is constrained in terms of other variables.
                          // We assume that all variables in this constraint have the same FE type (this
                          // is asserted below), and hence we can create the constraint row contribution
                          // by multiplying their_dof_value by the corresponding row of the transformation
                          // matrix.

                          const std::set<unsigned int> & variables = periodic->get_variables();
                          neigh_dof_indices_all_variables.resize(variables.size());
                          unsigned int index = 0;
                          for(unsigned int other_var : variables)
                            {
                              libmesh_assert_msg(base_fe_type == dof_map.variable_type(other_var), "FE types must match for all variables involved in constraint");

                              Real var_weighting = periodic->get_transformation_matrix()(variable_number, other_var);
                              constraint_row->emplace(neigh_dof_indices_all_variables[index][i],
                                                      var_weighting*their_dof_value);
                              index++;
                            }
                        }

                    }
                }
            }
          // p refinement constraints:
          // constrain dofs shared between
          // active elements and neighbors with
          // lower polynomial degrees
#ifdef LIBMESH_ENABLE_AMR
          const unsigned int min_p_level =
            neigh->min_p_level_by_neighbor(elem, elem->p_level());
          if (min_p_level < elem->p_level())
            {
              // Adaptive p refinement of non-hierarchic bases will
              // require more coding
              libmesh_assert(my_fe->is_hierarchic());
              dof_map.constrain_p_dofs(variable_number, elem,
                                       s, min_p_level);
            }
#endif // #ifdef LIBMESH_ENABLE_AMR
        }
    }
}

#endif // LIBMESH_ENABLE_PERIODIC

// ------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT FEGenericBase<Real>;
template class LIBMESH_EXPORT FEGenericBase<RealGradient>;

} // namespace libMesh
