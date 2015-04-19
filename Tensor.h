#ifndef TENSOR_H_
#define TENSOR_H_

#include <vector>

#include "Dimension.h"
#include "ContractionIterator.h"

template <class T>
struct Tensor {
  typedef std::vector<Dimension> Geometry;
  typedef std::vector<unsigned int> Coordinates;
  T * elements;
  Geometry geometry;
  unsigned int size();
  T & at(Coordinates &);
  T & at(ContractionIterator &);

  void zero();

  Tensor(Geometry &);

  static Tensor<T> contract(std::vector<Tensor<T> > &);
};

#include "Tensor.tcc"

#endif
