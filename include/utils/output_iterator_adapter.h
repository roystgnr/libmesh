// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_OUTPUT_ITERATOR_ADAPTER_H
#define LIBMESH_OUTPUT_ITERATOR_ADAPTER_H

// Local includes

// C++ includes
#include <iterator>

namespace libMesh
{

// A class for templated methods that expect packed_range compatible
// output iterator arguments but need to be handed standard C++ output
// iterators instead.
//
// The default version of the class inserts object values, relying on
// the inner iterator to copy them, then deletes the heap-allocated
// originals.
template <typename T, typename InnerIterator>
struct output_iterator_adapter
  : std::iterator<std::output_iterator_tag, T>
{
  output_iterator_adapter(InnerIterator i) : _iter(i) {}

  void operator=(T* obj) {
    *_iter = *obj;
    delete obj;
  }

  output_iterator_adapter& operator++() {
    ++_iter;
    return *this;
  }

  output_iterator_adapter operator++(int) {
    output_iterator_adapter result(*this);
    ++(*this);
    return result;
  }

  // We don't return a reference-to-T here because we don't want to
  // construct one or have any of its methods called.
  output_iterator_adapter& operator*() { return *this; }

private:
  InnerIterator _iter;
};


// The pointer-specialized version of the class inserts object
// pointers, leaving it up to the user to handle eventual deallocation
// of the heap-allocated objects.
//
// We ought to have a unique_ptr-specialized version too, but we can't
// rely on it since our UniquePtr still falls back on the
// container-unfriendly AutoPtr in some configurations.
template <typename T, typename InnerIterator>
struct output_iterator_adapter
  : std::iterator<std::output_iterator_tag, T*>
{
  output_iterator_adapter(InnerIterator i) : _iter(i) {}

  void operator=(T* obj) {
    *_iter = obj;
  }

  output_iterator_adapter& operator++() {
    ++_iter;
    return *this;
  }

  output_iterator_adapter operator++(int) {
    output_iterator_adapter result(*this);
    ++(*this);
    return result;
  }

  // We don't return a reference-to-T here because we don't want to
  // construct one or have any of its methods called.
  output_iterator_adapter& operator*() { return *this; }

private:
  InnerIterator _iter;
};

} // namespace libMesh


#endif // LIBMESH_OUTPUT_ITERATOR_ADAPTER_H
