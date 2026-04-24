#include <libmesh/libmesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class AllRBBTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the Mesh Extruder
   * with the optional object callback for setting custom subdomain IDs.
   * We pass a custom object for generating subdomains based on the old element
   * ID and the current layer and assert the proper values.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( AllRBBTest );

  CPPUNIT_TEST( testAllRBBNodeElem );
  CPPUNIT_TEST( testAllRBBEdge );
  CPPUNIT_TEST( testAllRBBEdge3 );

  // 2D tests
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testAllRBBTri );
  CPPUNIT_TEST( testAllRBBTri6 );
  CPPUNIT_TEST( testAllRBBQuad );
  CPPUNIT_TEST( testAllRBBQuad8 );
  CPPUNIT_TEST( testAllRBBQuad9 );
#endif

  // 3D tests
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testAllRBBTet );
  CPPUNIT_TEST( testAllRBBTet10 );
  CPPUNIT_TEST( testAllRBBHex );
  CPPUNIT_TEST( testAllRBBHex20 );
  CPPUNIT_TEST( testAllRBBHex27 );
  CPPUNIT_TEST( testAllRBBPrism6 );
  CPPUNIT_TEST( testAllRBBPrism18 );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:
  // We don't do anything interesting to affine elements in all_rbb(),
  // but we can verify that we're not screwing them up.
  void test_box(ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    auto dim = Elem::type_to_dim_map[elem_type];

    MeshTools::Generation::build_cube(mesh,
                                        2, dim > 1 ? 1 : 0, dim > 2 ? 1 : 0,
                                        0., 1.,
                                        0., 1.,
                                        0., 1.,
                                        elem_type);

    MeshTools::Modification::all_rbb(mesh);

    CPPUNIT_ASSERT_EQUAL(dof_id_type(2), mesh.n_elem());

    unsigned char weight_index = mesh.default_mapping_data();

    for (auto & elem : mesh.element_ptr_range())
    {
      CPPUNIT_ASSERT_EQUAL(elem->mapping_type(), RATIONAL_BERNSTEIN_MAP);
      CPPUNIT_ASSERT(elem->has_affine_map());
    }

    for (auto & node : mesh.node_ptr_range())
    {
      const Real w = node->get_extra_datum<Real>(weight_index);

      CPPUNIT_ASSERT_EQUAL(Real(1), w);
    }
  }

public:
  void setUp() {}

  void tearDown() {}

  void testAllRBBNodeElem() { LOG_UNIT_TEST; test_box(NODEELEM); }
  void testAllRBBEdge() { LOG_UNIT_TEST; test_box(EDGE2); }
  void testAllRBBEdge3() { LOG_UNIT_TEST; test_box(EDGE3); }
  void testAllRBBTri() { LOG_UNIT_TEST; test_box(TRI3); }
  void testAllRBBTri6() { LOG_UNIT_TEST; test_box(TRI6); }
  void testAllRBBQuad() { LOG_UNIT_TEST; test_box(QUAD4); }
  void testAllRBBQuad8() { LOG_UNIT_TEST; test_box(QUAD8); }
  void testAllRBBQuad9() { LOG_UNIT_TEST; test_box(QUAD9); }
  void testAllRBBTet() { LOG_UNIT_TEST; test_box(TET4); }
  void testAllRBBTet10() { LOG_UNIT_TEST; test_box(TET10); }
  void testAllRBBHex() { LOG_UNIT_TEST; test_box(HEX8); }
  void testAllRBBHex20() { LOG_UNIT_TEST; test_box(HEX20); }
  void testAllRBBHex27() { LOG_UNIT_TEST; test_box(HEX27); }
  void testAllRBBPrism6() { LOG_UNIT_TEST; test_box(PRISM6); }
  void testAllRBBPrism18() { LOG_UNIT_TEST; test_box(PRISM18); }

  // We still don't support general Polys, Tri7s or anything with them
  // as faces, infinite elements, or anything above quadratic.
};


CPPUNIT_TEST_SUITE_REGISTRATION( AllRBBTest );
