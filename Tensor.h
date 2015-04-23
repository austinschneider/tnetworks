#ifndef TENSOR_H_
#define TENSOR_H_

#include <vector>

#include "Dimension.h"
#include "ContractionIterator.h"

typedef std::vector<Dimension> Geometry;
typedef std::vector<unsigned int> Coordinates;

template <class T>
struct Tensor {
  T * elements;
  Geometry geometry;
  unsigned int size();
  T & at(Coordinates &);
  T & at(ContractionIterator &);

  void zero();

  Tensor();
  Tensor(Geometry &);
  Tensor(const Tensor<T> &);
  ~Tensor();

  static Tensor<T> contract(std::vector<Tensor<T> > &);
};

#include "Tensor.tcc"

#endif
